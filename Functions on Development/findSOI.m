function [rsoi] = findSOI(dist, mu1, mu2)

% DESCRIZIONE : 
% questa funzione permette di calcolare il raggio della sfera di influenza
% di un corpo secondario rispetto al corpo primario
%
% INPUT : 
% - dist : distanza tra il primario e il secondario (km)
% - mu1  : parametro gravitazionale del primario (km3/s2)
% - mu2  : parametro gravitazionale del secondario (km3/s2)
%
% OUTPUT :
% - rsoi : raggio della sfera di influenza del secondario (km)
%
% -------------------------------------------------------------------------


rsoi = dist*(mu2/mu1)^(2/5);

end