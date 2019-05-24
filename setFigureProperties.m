function f = setFigureProperties

% This function aims to set figure properties

% INPUT: -
% OUTPUT: f – figure handle

   % prepare environment for plot
   f = figure('Visible','Off');
   set(f,'PaperPositionMode','auto');
   set(f,'units','normalized','outerposition',[0 0 1 1]);
   set(f, 'PaperType', 'a4');
   set(f, 'PaperOrientation', 'Portrait');
end