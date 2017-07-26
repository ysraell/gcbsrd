function res=normalizar(entrada)
    sinal=entrada;
    res=(sinal-min(min(sinal)))/(max(max(sinal))-min(min(sinal)));
end