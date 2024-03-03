#%%
import os

MAIN_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_PATH = os.path.join(MAIN_PATH, 'data')
CLIMA_PATH = os.path.join(DATA_PATH, '01_raw')
DB_PATH = os.path.join(DATA_PATH, '02_intermediate')

NOTE_PATH = os.path.join(MAIN_PATH, 'notebooks')
CODE_PATH = os.path.join(MAIN_PATH, 'src')


#%%
# Função para ler os dados do CSV e obter a região
def obter_regiao(arquivo_csv):
    with open(arquivo_csv, 'r', encoding='latin-1') as f:
        # Lendo a primeira linha para obter a região
        primeira_linha = f.readline().strip()
        # Dividindo a linha pelo separador ":;"
        regiao = primeira_linha.split(":;")[-1]
    return regiao

# Função para ler os dados do CSV e obter o estado
def obter_estado(arquivo_csv):
    with open(arquivo_csv, 'r', encoding='latin-1') as f:
        # Lendo todas as linhas do arquivo
        linhas = f.readlines()
        # Pegando a segunda linha para obter o estado
        segunda_linha = linhas[1].strip()
        # Dividindo a linha pelo separador ":;"
        estado = segunda_linha.split(":;")[-1]
    return estado

#%%
# Diretório raiz onde estão as pastas por ano
diretorio_raiz = CLIMA_PATH

# Listando todos os anos
anos = [str(ano) for ano in range(2000, 2025)]

#%%
# Para cada ano
for ano in anos:
    # Caminho para a pasta do ano
    pasta_ano = os.path.join(diretorio_raiz, ano)
    # Listando todos os arquivos CSV na pasta do ano
    arquivos_csv = [arquivo for arquivo in os.listdir(pasta_ano) if arquivo.endswith('.CSV')]

    # Para cada arquivo CSV
    for arquivo_csv in arquivos_csv:
        
        # Obtendo a região do arquivo
        regiao = obter_regiao(os.path.join(pasta_ano, arquivo_csv))
        
        # Criando o diretório da região se não existir
        pasta_regiao = os.path.join(diretorio_raiz, regiao)
        if not os.path.exists(pasta_regiao):
            os.makedirs(pasta_regiao)
        
        # Movendo o arquivo para o diretório da região
        novo_caminho = os.path.join(pasta_regiao, arquivo_csv)
        os.rename(os.path.join(pasta_ano, arquivo_csv), novo_caminho)

'''# Listando todas as pastas de região
pastas_regiao = [pasta for pasta in os.listdir(diretorio_raiz) if os.path.isdir(os.path.join(diretorio_raiz, pasta))]

# Para cada pasta de região
for pasta_regiao in pastas_regiao:
    # Caminho para a pasta da região
    pasta_regiao_atual = os.path.join(diretorio_raiz, pasta_regiao)
    
    # Listando todos os arquivos CSV na pasta da região
    arquivos_csv = [arquivo for arquivo in os.listdir(pasta_regiao_atual) if arquivo.endswith('.CSV')]
    
    # Para cada arquivo CSV
    for arquivo_csv in arquivos_csv:
        # Obtendo o estado do arquivo
        estado = obter_estado(os.path.join(pasta_regiao_atual, arquivo_csv))
        
        # Criando o diretório do estado dentro da pasta da região se não existir
        pasta_estado = os.path.join(pasta_regiao_atual, estado)
        if not os.path.exists(pasta_estado):
            os.makedirs(pasta_estado)
        
        # Movendo o arquivo para o diretório do estado
        novo_caminho = os.path.join(pasta_estado, arquivo_csv)
        os.rename(os.path.join(pasta_regiao_atual, arquivo_csv), novo_caminho)'''
