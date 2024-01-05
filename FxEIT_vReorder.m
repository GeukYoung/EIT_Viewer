function [idx_change] = FxEIT_vReorder(num_elec, mask)
    nProj = length(num_elec);
    if nargin < 2
        mask = FxEIT_mask(nProj);
    end
    if mask == 'x'
        mask = false(nProj,1);
    end
    
    idx_elec.origin = 1:nProj;
    idx_elec.change = num_elec;
    
    %% proj table
    for i = 1:nProj
        for j = 1:nProj
            inj1 = i-1; inj2 = i; mea1 = j-1; mea2 = j;
            if inj1 == 0
                inj1 = nProj;
            end
            if mea1 == 0
                mea1 = nProj;
            end
            table_proj.origin((i-1)*nProj+j,:) = idx_elec.origin([inj1 inj2 mea1 mea2]);
            table_proj.change((i-1)*nProj+j,:) = idx_elec.change([inj1 inj2 mea1 mea2]);
        end
    end
    % sort
    table_proj.origin = [sort(table_proj.origin(:,1:2),2) sort(table_proj.origin(:,3:4),2)];
    table_proj.change = [sort(table_proj.change(:,1:2),2) sort(table_proj.change(:,3:4),2)];
    % mask
    table_proj.origin(mask,:) = [];
    table_proj.change(mask,:) = [];
    
    %% find index (origin projnum) 
    idx_change = zeros(length(table_proj.origin),1);
    for i = 1:size(table_proj.origin,1)
        idx_change(i) = find(sum(abs(table_proj.origin(:,:) - table_proj.change(i,:)),2)==0);
    end

    if length(unique(idx_change)) ~= size(table_proj.origin,1)
        error('matching miss');
    end
    
%     figure; plot(idx_change,'.');
end

function [mask] = FxEIT_mask(ch,first_sat)
temp = zeros(ch,ch);
if nargin == 1
    first_sat = [ch,1,2];
end

if length(first_sat) == 4
    cnt1 = first_sat(1); cnt2 = first_sat(2); cnt3 = first_sat(3); cnt4 = first_sat(4);
    for i = 1:ch
        temp([cnt1,cnt2,cnt3,cnt4],i) = 1;
        cnt1 = cnt1 + 1;
        cnt2 = cnt2 + 1;
        cnt3 = cnt3 + 1;
        cnt4 = cnt4 + 1;
        if cnt1 > ch
            cnt1 = 1;
        end
        if cnt2 > ch
            cnt2 = 1;
        end
        if cnt3 > ch
            cnt3 = 1;
        end
        if cnt4 > ch
            cnt4 = 1;
        end
    end
else
    cnt1 = first_sat(1); cnt2 = first_sat(2); cnt3 = first_sat(3);
    for i = 1:ch
        temp([cnt1,cnt2,cnt3],i) = 1;
        cnt1 = cnt1 + 1;
        cnt2 = cnt2 + 1;
        cnt3 = cnt3 + 1;
        if cnt1 > ch
            cnt1 = 1;
        end
        if cnt2 > ch
            cnt2 = 1;
        end
        if cnt3 > ch
            cnt3 = 1;
        end
    end
end

temp2 = reshape(temp,ch*ch,1);
[mask, ~] = find(temp2 == 1);
end