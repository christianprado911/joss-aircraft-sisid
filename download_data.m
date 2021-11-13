%% Faz o download e descompacta os dados de voo
websave('flt_data.zip', ...
    'https://arc.aiaa.org/doi/suppl/10.2514/4.102790/suppl_file/flt_data.zip')
unzip('flt_data.zip');
