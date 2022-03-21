function stop = WriteToFile1(t,y,Flag)
stop = false;
persistent fileout

if strcmp(Flag,'init')
  FileName = 'VdP_Result_1.txt';
  fileout = fopen(FileName,'w');  
elseif strcmp(Flag,'') 
  fprintf(fileout,'  %18.16e',t);
  for l = 1:size(y,2)
    fprintf(fileout,'  %18.16e',y(l));
  end 
  fprintf(fileout,'\r\n');
elseif strcmp(Flag,'done')
  fclose(fileout);
end

