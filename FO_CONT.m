function [prob]=FO_CONT(prob,TO_ref)
        
    % create row of max ToF
    prob=EO_t0CONT(prob,TO_ref);
    
    % continuate along each column
    L=length(prob);

    for id=1:L

        prob=EO_tfCONT(prob,TO_ref,id);

    end

end