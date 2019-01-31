function y=changeFs(x,Fini,Fend)


    [N,D]=rat(Fend/Fini);
    
    xup=interp(x,N);
  
    y=xup(1:D:end);