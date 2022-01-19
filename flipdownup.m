function Y = FlipDownUp(N,num,i,j)
 
 %Y = bitxor(num, 2^i+2^j);
 Y=0;
 
 for ii = 0 : N-1
     
    if ii==i
      Y = Y + 0*3^ii;
    elseif ii==j 
      Y = Y + 2*3^ii;  
    else
      Y = Y + tritget(num,ii+1)*3^ii;
    end  
 
end
