function [ sim ] = combineSim3( sim1, sim2 ,a)
%COMBINESIM Summary of this function goes here
%   Detailed explanation goes here

[rows, cols]=size(sim1);

sim=zeros(rows, cols);

 for i=1:rows
     for j=1:cols
        
         
         if (sim1(i,j)==0)
             sim(i,j)=sim2(i,j);
        else
             %sim(i,j)=sim1(i,j);
             sim(i,j)=(a*sim1(i,j)+(1-a)*sim2(i,j));
         end
       
                
                     
     end
 end

end

