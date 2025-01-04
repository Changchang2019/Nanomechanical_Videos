% decimate properly
function [x_decim] = decim_joel(x,sampsPerSym,decim_offset)
x_decim = nan(1,floor(size(x,1)/sampsPerSym));
   for n=1:floor(size(x,1)/sampsPerSym)
           if 1+(n-1)*sampsPerSym + decim_offset<= size(x,1)
                x_decim(n) = x(1+(n-1)*sampsPerSym + decim_offset);
           else
               x_decim(n) = x(1+(n-1)*sampsPerSym + 0);
           end
   end
end