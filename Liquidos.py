from Constants import liquidos

for liquido in liquidos:
    dados = pd.read_excel(liquido)
    dados = dados.sort_values('temp [K]')
    x = np.array(dados['temp [K]'])
    y_log = np.array(dados['visc [log10 Pa.s]'])