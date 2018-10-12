function [rindex] = kkre_hilbert(absorbance )
    rindex = -imag( hilbert( absorbance ) ); 
end

