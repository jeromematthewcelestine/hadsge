classdef MCMCResults < handle
properties
    n_chains
    n_rows_per_chain
    n_rows_per_raw_chain
    acceptance_rate
    n_rows_total
    
    n_params
    
    chains_raw
    chains
    lengths_raw
    
    combined
    
    Rhat
end
methods
    function o = MCMCResults(path, fraction_to_keep)
        file_list = dir([path,'/chain*']);
        o.n_chains = length(file_list);

        o.chains_raw = cell(o.n_chains,1);
        
        o.combined = [];

        acceptances = 0;
        o.n_rows_total = 0;

        for chain_idx = 1:o.n_chains
            table_name = file_list(chain_idx).name;
            raw_chain = dlmread([path,'/',table_name]);
            o.chains_raw{chain_idx} = raw_chain;
            o.lengths_raw(chain_idx) = size(raw_chain,1);
        end
        o.n_rows_per_raw_chain = min(o.lengths_raw');
        
        for chain_idx = 1:o.n_chains
            raw_chain = o.chains_raw{chain_idx};
            
            o.n_params = size(raw_chain,2) - 2;
            o.chains_raw{chain_idx} = raw_chain;
            
            o.n_rows_per_chain = floor(fraction_to_keep * o.n_rows_per_raw_chain);
            
            chain = raw_chain(1+(end - o.n_rows_per_chain):end,:);
            o.chains{chain_idx} = chain;
            
            if (isempty(o.combined))
                o.combined = chain;
            else
                o.combined = [o.combined; chain];
            end
            
            o.n_rows_per_chain = size(chain,1);
            o.n_rows_total = o.n_rows_total + o.n_rows_per_chain;

            acceptances = acceptances + sum(chain(:,end));
        end
        
        o.acceptance_rate = acceptances / o.n_rows_total;
    end
end
end