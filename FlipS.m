function Y = FlipS(N,num,i,j)
 
 %Y = bitxor(num, 2^i+2^j);
 Y=0;
 
 for ii = 0 : N-1
     
    if ii==i
      Y = Y + tritget(num,j+1)*3^ii;
    elseif ii==j 
      Y = Y + tritget(num,i+1)*3^ii;  
    else
      Y = Y + tritget(num,ii+1)*3^ii;
    end  
 
end
