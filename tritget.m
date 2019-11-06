function j = tritget(num,i)
 
 div = num;

 if (i>1) && (num <= 3^(i-1)-1)
    j=0;  
 else    
    for ii = 1 : i
       j = rem(div,3);
       div = round( (div - j)/3 );
    end   
 end
    
end
