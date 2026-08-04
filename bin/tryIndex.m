function val = tryIndex(mat, iLon1, iLat1, iLat2, iLon2)
    % If mat is 2D and size matches [length(lons), length(lats)] assume (lon,lat)
    sz = size(mat);
    val = NaN;
    if numel(sz) == 2
        if sz(1) >= iLon1 && sz(2) >= iLat1
            try val = mat(iLon1, iLat1); return; end
        end
        if sz(1) >= iLat2 && sz(2) >= iLon2
            try val = mat(iLat2, iLon2); return; end
        end
        % fallback: linear index if sizes match total
        if numel(mat) >= max(iLon1,iLat1)
            try val = mat(iLon1); return; end
        end
    elseif numel(sz) == 3
        % if there's an extra singleton or layer dim, pick the first page
        % try common permutations
        try
            if sz(1) >= iLon1 && sz(2) >= iLat1
                val = mat(iLon1, iLat1, 1); return;
            elseif sz(1) >= iLat2 && sz(2) >= iLon2
                val = mat(iLat2, iLon2, 1); return;
            end
        catch
        end
    end
    % if all fails, return NaN
end