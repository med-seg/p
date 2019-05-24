function printFigureWithName(f,filename)
   pictureFileName = strrep(strrep(filename,sprintf('\r\n'),'_'),sprintf('\n'),'_');
   print('-dpng','-r0',pictureFileName); % save to file
   close(f);
end