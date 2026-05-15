'''
file.py — Módulo de Leitura de Arquivos de Micro-FTIR
======================================================
Leitura nativa dos principais formatos de imagens
hiperespectrais geradas por microscópios de FTIR.

FORMATOS SUPORTADOS:
  fsm(arq)          — Perkin Elmer Spotlight (.fsm)  [parser binário nativo]
  age(path)         — Agilent (.dmd + .dmt)
  npz_save/npz_load — formato interno numpy (.npz)

─────────────────────────────────────────────────────
ESTRUTURA INTERNA DO DICIONÁRIO HSP (padrão da lib)
─────────────────────────────────────────────────────
Todos os leitores retornam um dicionário com as chaves:
    r        — ndarray float32 (n_espectros × n_pontos): absorbâncias
    wn       — ndarray float64 (n_pontos,): números de onda em cm⁻¹
    dx       — int: número de linhas da imagem (eixo vertical)
    dy       — int: número de colunas da imagem (eixo horizontal)
    sel      — ndarray bool (n_espectros,): máscara de pixels válidos
    filename — str: nome do arquivo fonte
    log      — str: histórico de operações aplicadas aos dados

─────────────────────────────────────────────────────
ESTRUTURA BINÁRIA DO FSM (Perkin Elmer Spotlight)
─────────────────────────────────────────────────────
  bytes 0–3   : assinatura 'PEPE' (identificador do formato)
  bytes 4–43  : descrição do experimento (40 bytes, codificação UTF-8)
  bytes 44+   : sequência de blocos com estrutura:
                  [2 bytes] block_id   (unsigned short, little-endian)
                  [4 bytes] block_size (signed int,    little-endian)
                  [N bytes] dados do bloco

  Blocos relevantes:
    ID 5100 : metadados do mapeamento (n_x, n_y, n_z, z_start, z_end, z_delta)
    ID 5104 : informações do instrumento (analista, data, modelo, etc.)
    ID 5105 : dados espectrais brutos (array float32 de todos os espectros)

  Convenção de eixos do FSM:
    n_x → número de colunas da imagem → dy na convenção da biblioteca
    n_y → número de linhas  da imagem → dx na convenção da biblioteca
'''

import struct
import os
import glob
import numpy as np


# ─────────────────────────────────────────────────────────────────────────────
# FUNÇÕES INTERNAS DE PARSING BINÁRIO DO FSM
# ─────────────────────────────────────────────────────────────────────────────

def _fsm_block_info(data):
    '''
    Lê o cabeçalho de 6 bytes de um bloco FSM.

    O cabeçalho é composto por:
        2 bytes — block_id   (unsigned short, little-endian): tipo do bloco
        4 bytes — block_size (signed int,    little-endian): tamanho em bytes

    Parâmetros
    ----------
    data : bytes
        Exatamente 6 bytes lidos da posição atual do arquivo.

    Retorna
    -------
    (block_id, block_size) : tuple(int, int)
    '''
    return struct.unpack('<Hi', data)   # '<' = little-endian, 'H' = uint16, 'i' = int32


def _fsm_decode_5100(data):
    '''
    Decodifica o bloco FSM de ID 5100: metadados principais do mapeamento.

    Este bloco contém todas as informações geométricas e espectrais necessárias
    para reconstruir a imagem hiperespectral:
        - Dimensões espaciais: n_x (colunas), n_y (linhas)
        - Número de pontos espectrais: n_z
        - Intervalo espectral: z_start (cm⁻¹ inicial), z_end (cm⁻¹ final),
          z_delta (passo entre pontos)
        - Espaçamento espacial: x_delta, y_delta (µm por pixel)
        - Posição inicial do mapeamento: x_init, y_init

    Estrutura do bloco após o campo de nome:
        104 bytes com o formato '<ddddddddddiiihBhBhBhB'
        (combinação de doubles, ints e bytes, little-endian)

    Parâmetros
    ----------
    data : bytes
        Conteúdo bruto do bloco 5100.

    Retorna
    -------
    dict com campos: name, x_delta, y_delta, z_delta, z_start, z_end,
                     x_init, y_init, n_x, n_y, n_z, resolution, transmission
    '''
    # Os primeiros 2 bytes informam o tamanho do nome do experimento
    name_size = struct.unpack('<h', data[:2])[0]
    name = data[2:name_size + 2].decode('utf8', errors='replace')

    # Formato dos 104 bytes de metadados que seguem o nome
    fmt = '<ddddddddddiiihBhBhBhB'
    header_size = 104
    (x_delta, y_delta, z_delta, z_start, z_end,
     z_4d_start, z_4d_end, x_init, y_init, z_init,
     n_x, n_y, n_z,
     _, text1, _, text2, resolution, text3, transmission, text4
     ) = struct.unpack(fmt, data[name_size + 2: name_size + header_size + 2])

    return {
        'name':         name,
        'x_delta':      x_delta,      # µm por pixel (eixo x)
        'y_delta':      y_delta,      # µm por pixel (eixo y)
        'z_delta':      z_delta,      # passo espectral em cm⁻¹
        'z_start':      z_start,      # número de onda inicial em cm⁻¹
        'z_end':        z_end,        # número de onda final em cm⁻¹
        'x_init':       x_init,       # posição x inicial do mapeamento
        'y_init':       y_init,       # posição y inicial do mapeamento
        'n_x':          n_x,          # número de colunas
        'n_y':          n_y,          # número de linhas
        'n_z':          n_z,          # número de pontos espectrais
        'resolution':   resolution,   # resolução espectral em cm⁻¹
        'transmission': transmission, # modo de medição (transmissão ou reflexão)
    }


def _fsm_decode_5104(data):
    '''
    Decodifica o bloco FSM de ID 5104: informações do instrumento.

    Este bloco contém texto estruturado com metadados operacionais como:
        - Nome do analista
        - Data/hora da medição
        - Nome da imagem
        - Modelo do instrumento
        - Número de série
        - Versão do software

    O bloco é codificado como sequência de registros com tags de 2 bytes:
        '#u' — string: 2 bytes de tamanho + N bytes UTF-8 + 6 bytes de padding
        '$u' — inteiro curto: 2 bytes + 6 bytes de padding
        ',u' — inteiro curto: 2 bytes (sem padding)

    Parâmetros
    ----------
    data : bytes
        Conteúdo bruto do bloco 5104.

    Retorna
    -------
    dict com campos disponíveis entre: analyst, date, image_name,
    instrument_model, instrument_serial_number, instrument_software_version
    '''
    text = []
    i = 0
    while i + 2 < len(data):
        tag = data[i:i + 2]
        if tag == b'#u':
            # String UTF-8: lê 2 bytes de tamanho, depois N bytes de conteúdo
            i += 2
            size = struct.unpack('<h', data[i:i + 2])[0]
            i += 2
            text.append(data[i:i + size].decode('utf8', errors='replace'))
            i += size + 6   # pula 6 bytes de padding após a string
        elif tag in (b'$u', b',u'):
            # Inteiro curto (2 bytes)
            i += 2
            text.append(struct.unpack('<h', data[i:i + 2])[0])
            i += 2
            if tag == b'$u':
                i += 6   # padding apenas para '$u'
        else:
            i += 1   # byte desconhecido: avança um byte

    # Mapeia os valores extraídos em ordem para os campos conhecidos
    info = {}
    labels = ['analyst', None, 'date', None, 'image_name',
              'instrument_model', 'instrument_serial_number',
              'instrument_software_version']
    for idx, label in enumerate(labels):
        if label and idx < len(text):
            info[label] = text[idx]
    return info


def _fsm_decode_5105(data):
    '''
    Decodifica o bloco FSM de ID 5105: dados espectrais brutos.

    Este bloco contém todos os espectros do mapeamento concatenados
    sequencialmente como valores float32 em precisão simples. A ordem
    segue o padrão de varredura do Spotlight: linha por linha, da esquerda
    para a direita.

    Parâmetros
    ----------
    data : bytes
        Conteúdo bruto do bloco 5105.

    Retorna
    -------
    ndarray float32 : vetor 1D com todos os valores espectrais brutos.
        Reshape posterior para (n_espectros, n_pontos) é feito em _read_fsm_binary.
    '''
    return np.frombuffer(data, dtype=np.float32)


# Dicionário de decoders mapeando block_id → função de decodificação
_FSM_DECODERS = {
    5100: _fsm_decode_5100,   # metadados geométricos e espectrais
    5104: _fsm_decode_5104,   # informações do instrumento
    5105: _fsm_decode_5105,   # dados espectrais brutos
}


def _read_fsm_binary(path):
    '''
    Função principal de parsing do arquivo FSM binário.

    Lê o arquivo completo na memória e percorre todos os blocos sequencialmente,
    decodificando os relevantes (IDs 5100, 5104, 5105) e ignorando os demais.

    Parâmetros
    ----------
    path : str
        Caminho completo para o arquivo .fsm.

    Retorna
    -------
    amplitudes : ndarray float32, shape (n_espectros, n_pontos)
        Dados espectrais brutos em unidades de transmitância (%).
    wavelength : ndarray float64, shape (n_pontos,)
        Vetor de números de onda em cm⁻¹, calculado a partir de z_start,
        z_end e z_delta lidos do bloco 5100.
    meta : dict
        Dicionário com todos os metadados extraídos dos blocos 5100 e 5104,
        mais os campos 'signature', 'description' e 'filename'.
    '''
    with open(path, 'rb') as f:
        content = f.read()   # carrega o arquivo inteiro na memória

    # Valida a assinatura do formato (primeiros 4 bytes devem ser 'PEPE')
    if content[:4] != b'PEPE':
        raise ValueError(
            'Arquivo não é um FSM válido — assinatura esperada: PEPE\n'
            'Arquivo: ' + path
        )

    # Metadados iniciais lidos do cabeçalho fixo
    meta = {
        'signature':   content[:4],
        'description': content[4:44].decode('utf8', errors='replace').strip('\x00'),
        'filename':    os.path.basename(path),
    }

    spectral_blocks = []   # acumula os arrays float32 do bloco 5105
    pos = 44               # posição inicial: após os 44 bytes de cabeçalho

    # Percorre todos os blocos do arquivo
    while pos + 6 < len(content):
        # Lê cabeçalho do bloco: ID (2 bytes) + tamanho (4 bytes)
        block_id, block_size = _fsm_block_info(content[pos:pos + 6])
        pos += 6
        block_data = content[pos:pos + block_size]
        pos += block_size

        # Decodifica o bloco se existir um decoder para seu ID
        if block_id in _FSM_DECODERS:
            decoded = _FSM_DECODERS[block_id](block_data)
            if isinstance(decoded, dict):
                # Blocos 5100 e 5104: atualiza metadados
                meta.update(decoded)
            else:
                # Bloco 5105: acumula array de espectros brutos
                spectral_blocks.append(decoded)

    # Empilha todos os blocos de espectros em um único array 2D
    amplitudes = np.squeeze(np.array(spectral_blocks, dtype=np.float32))

    # Reconstrói o vetor de números de onda a partir dos metadados do bloco 5100
    wavelength = np.arange(
        meta['z_start'],
        meta['z_end'] + meta['z_delta'],
        meta['z_delta'],
        dtype=np.float64
    )
    # Garante que o vetor tenha exatamente n_z pontos (evita erros de arredondamento)
    if 'n_z' in meta:
        wavelength = wavelength[:meta['n_z']]

    return amplitudes, wavelength, meta


# ─────────────────────────────────────────────────────────────────────────────
# API PÚBLICA — LEITURA FSM
# ─────────────────────────────────────────────────────────────────────────────

def fsm(arq):
    '''
    Lê um arquivo FSM da Perkin Elmer Spotlight e retorna o dicionário hsp.

    Converte automaticamente os dados de transmitância (%) para absorbância
    usando a Lei de Beer-Lambert:
        A = -log₁₀(T/100)

    O vetor de números de onda é invertido para seguir a convenção padrão
    de espectroscopia (valores decrescentes da esquerda para a direita).

    Parâmetros
    ----------
    arq : str
        Caminho completo para o arquivo .fsm a ser lido.

    Retorna
    -------
    data : dict hsp com as chaves:
        r        — ndarray float32 (n_espectros × n_pontos): absorbâncias
        wn       — ndarray float64 (n_pontos,): números de onda em cm⁻¹
        dx       — int: número de linhas da imagem (eixo vertical)
        dy       — int: número de colunas da imagem (eixo horizontal)
        sel      — ndarray bool (n_espectros,): todos True na leitura inicial
        filename — str: nome do arquivo .fsm
        log      — str: mensagem de log inicial
    '''
    amplitudes, wavelength, meta = _read_fsm_binary(arq)

    # Convenção de eixos do FSM:
    #   n_x = número de colunas → dy na convenção da biblioteca (eixo horizontal)
    #   n_y = número de linhas  → dx na convenção da biblioteca (eixo vertical)
    n_x = int(meta.get('n_x', 1))
    n_y = int(meta.get('n_y', 1))
    dx  = n_y   # linhas da imagem
    dy  = n_x   # colunas da imagem

    # Garante shape 2D mesmo quando há apenas 1 espectro (ex.: aquisição pontual)
    if amplitudes.ndim == 1:
        amplitudes = amplitudes.reshape(1, -1)

    # Converte transmitância (%) → absorbância: A = -log₁₀(T/100)
    # np.fliplr inverte a ordem espectral (FSM armazena do maior para menor cm⁻¹)
    # np.clip evita log(0) ou log(negativo) em pixels ruidosos (limita mínimo a 1e-10)
    r = np.fliplr(amplitudes.astype('float32'))
    r = -np.log10(np.clip(0.01 * r, 1e-10, None))

    # Inverte o vetor de números de onda para corresponder à inversão de r
    wn = np.flipud(wavelength).astype('float64')

    data = {
        'r':        r,
        'wn':       wn,
        'dx':       np.array(dx),
        'dy':       np.array(dy),
        'filename': np.array(meta.get('filename', os.path.basename(arq))),
        'sel':      np.ones(r.shape[0], dtype=bool),   # todos os pixels válidos
        'log':      np.array('abrindo arquivo FSM: ' + os.path.basename(arq)),
    }

    print(' abrindo o arquivo  ' + arq)
    print('  >> ' + str(r.shape[0]) + ' espectros | ' +
          str(r.shape[1]) + ' pontos | ' +
          str(dx) + ' linhas x ' + str(dy) + ' colunas')
    return data


def get_fsm_files(path):
    '''
    Retorna uma lista de caminhos de todos os arquivos .fsm em um diretório.

    Muda o diretório de trabalho atual para path (efeito colateral via os.chdir).

    Parâmetros
    ----------
    path : str
        Caminho do diretório a ser pesquisado.

    Retorna
    -------
    list of str
        Lista com os caminhos completos de todos os arquivos .fsm encontrados.
    '''
    os.chdir(path)
    return [os.path.join(path, f) for f in glob.glob('*.fsm')]


# ─────────────────────────────────────────────────────────────────────────────
# API PÚBLICA — LEITURA AGILENT (.dmd + .dmt)
# ─────────────────────────────────────────────────────────────────────────────

def age(path):
    '''
    Lê um mosaico hiperespectral Agilent a partir de arquivos .dmd e .dmt.

    O microscópio Agilent divide a imagem em tiles de 32×32 pixels, cada um
    salvo como um arquivo .dmd separado. O arquivo .dmt contém os metadados
    espectrais (número de pontos, número de onda inicial e passo).

    Convenção de nomenclatura dos arquivos .dmd:
        <path>XXXX_YYYY.dmd
        onde XXXX = índice de tile no eixo x (coluna de tiles)
              YYYY = índice de tile no eixo y (linha de tiles)

    Estrutura do .dmt (binário):
        Offset  559 (int32)  : n_z (número de pontos espectrais)
        Offset  557 (int32)  : z_start (número de onda inicial × 1/steps)
        Offset  389 (float64): steps (passo espectral em cm⁻¹)

    Estrutura de cada .dmd (binário):
        255 floats de cabeçalho (ignorados)
        Restante: n_z × 32 × 32 floats → reorganizados para (32, 32, n_z)

    Parâmetros
    ----------
    path : str
        Caminho base do mosaico (sem o nome dos arquivos individuais).
        Ex: '/dados/experimento_01/' (os .dmd e .dmt devem estar nessa pasta)

    Retorna
    -------
    dados : dict hsp com as chaves:
        r        — ndarray float32 (n_espectros × n_pontos): absorbâncias
        wn       — ndarray float64 (n_pontos,): números de onda em cm⁻¹
        dx       — int: número de linhas totais (32 × n_tiles_x)
        dy       — int: número de colunas totais (32 × n_tiles_y)
        sel      — ndarray bool (n_espectros,): todos True na leitura inicial
        filename — str: path base do mosaico
        log      — str: mensagem de log inicial
    '''
    # Determina as dimensões do mosaico a partir dos nomes dos arquivos .dmd
    dx_list, dy_list = [], []
    for fname in glob.glob(path + '*.dmd'):
        # Extrai índices XXXX e YYYY do final do nome do arquivo
        parts = fname[-13:-4].split('_')
        dx_list.append(int(parts[0]))
        dy_list.append(int(parts[1]))

    # Número de tiles em cada eixo (+1 porque os índices começam em 0)
    dx = np.array(dx_list).max() + 1
    dy = np.array(dy_list).max() + 1

    dados = {
        'dx':       np.array(32 * dx),    # pixels totais no eixo x
        'dy':       np.array(32 * dy),    # pixels totais no eixo y
        'filename': path,
    }

    # Lê metadados espectrais do arquivo .dmt
    arq_dmt = glob.glob(path + '*.dmt')[0]
    tmp_i4  = np.fromfile(arq_dmt, dtype='i4')
    npoints = tmp_i4[559]       # número de pontos espectrais
    start   = tmp_i4[557]       # índice do número de onda inicial
    tmp_f64 = np.fromfile(arq_dmt)
    steps   = tmp_f64[389]      # passo espectral em cm⁻¹

    # Reconstrói o vetor de números de onda
    dados['wn'] = np.linspace(start * steps,
                               (start + npoints) * steps,
                               npoints)

    # Aloca array 3D para todo o mosaico: (linhas_totais, colunas_totais, n_pontos)
    r = np.zeros((int(32 * dx), int(32 * dy), int(npoints)), dtype='float32')

    # Lê e posiciona cada tile no array global
    for fname in glob.glob(path + '*.dmd'):
        y    = int(fname[-8:-4])    # índice do tile na direção y (colunas)
        x    = int(fname[-13:-9])   # índice do tile na direção x (linhas)
        data = np.fromfile(fname, 'f4')[255:]   # pula os 255 floats de cabeçalho
        # Reorganiza de (n_z × 32 × 32) para (32, 32, n_z)
        data = np.reshape(data, (-1, 32, 32))
        data = np.transpose(data, (2, 1, 0))
        data = data[:, ::-1, :]    # inverte eixo y para corrigir orientação
        # Insere o tile na posição correta do mosaico
        r[32 * x:32 * (x + 1), 32 * y:32 * (y + 1), :] = data

    # Converte para 2D (n_espectros × n_pontos) e cria máscara de pixels válidos
    dados['r']   = r.reshape(1024 * dx * dy, -1)
    dados['sel'] = np.ones(dados['r'].shape[0], dtype=bool)
    dados['log'] = np.array('abrindo mosaico agilent: ' + path)
    return dados


# ─────────────────────────────────────────────────────────────────────────────
# API PÚBLICA — SALVAR/CARREGAR FORMATO INTERNO (.npz)
# ─────────────────────────────────────────────────────────────────────────────

def npz_save(arq, data):
    '''
    Salva o dicionário hsp em formato .npz (arquivo comprimido NumPy).

    O formato .npz permite salvar e recuperar todo o estado de processamento
    (espectros, metadados, máscaras, log) em um único arquivo binário,
    sem precisar reprocessar o arquivo original.

    Parâmetros
    ----------
    arq : str
        Caminho de destino do arquivo (sem extensão, que será .npz).
    data : dict hsp
        Dicionário com os dados a salvar. Todas as chaves (r, wn, dx, dy,
        sel, filename, log, etc.) são preservadas.
    '''
    np.savez(arq, **data)


def npz_load(arq):
    '''
    Carrega um dicionário hsp previamente salvo em formato .npz.

    Parâmetros
    ----------
    arq : str
        Caminho do arquivo .npz a carregar.

    Retorna
    -------
    dict
        Dicionário com as mesmas chaves e arrays salvos por npz_save().
        allow_pickle=True é necessário para carregar strings e arrays
        de objetos (como 'filename' e 'log').
    '''
    f = np.load(arq, allow_pickle=True)
    return {k: f[k] for k in f.files}
