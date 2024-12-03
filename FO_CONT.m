function [prob]=FO_CONT(prob,TO_ref)
        
    % create row of max ToF
    prob=EO_t0CONT(prob,TO_ref);
    
    % continuate along each column
    L=length(prob);

    for id=1:L

        prob=EO_tfCONT(prob,TO_ref,id);
        prob(end).isTO=1;
        prob(end)=E2F_CONT(prob(end),id,L);

    end

    % for i=1:length(prob)
    %     prob(i)=E2F_CONT(prob(i),i,length(prob));
    % end

end

% t0 tf and e2f all only used in fo_cnt