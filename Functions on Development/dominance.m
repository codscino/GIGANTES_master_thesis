function dom=dominance(f,h);
% 
%   dominance
%   
%  [dom]=dominance(f,h)
%
%  Evaluates the dominance of  element i
%  over any other element j in f
%
%   f(i,:)>f(j,:) 
%
%   if h=0
%      simple dominance is computed, i.e. 
%      i is dominated by j if all the not equal components of j are better than i
%   if h=1
%      strong dominance, i.e.
%      i is dominante by j if any component of j is better than i and the
%      dominance level is the sum of all the better components
%
%   (c)  Massimiliano Vasile  2003
%
[nf,mf]=size(f);

dom=zeros(1,nf);    

if (mf==1)
    dom=f';
else
    if h<2
    for i=1:nf  
        for j=1:nf       
            if (i~=j)
                if (h==0)
                    nz=find(f(i,:)~=f(j,:));
                    if (all(f(i,nz)>f(j,nz)))&&(isempty(nz)==0)            
                        dom(i)=dom(i)+1;    
                    end      
                elseif (h==1)
                    if (any(f(i,:)>f(j,:)))
                        dom(i)=dom(i)+sum(f(i,:)>f(j,:));
                    end
                end
            end
        end        
    end
    else
        dom=sum([f(1,:)==f(2,:)])*max([f(1,:)>f(2,:)])+sum(f(1,:)>f(2,:));
    end
end

return