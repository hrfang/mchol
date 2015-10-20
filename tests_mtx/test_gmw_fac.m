fid = fopen('mtx_list.txt');

aline = fgetl(fid);
while ischar(aline)
    fileA = strtrim(aline);  % trim leading and tailing spaces etc.
    if length(fileA)
        fprintf('# factorizing %s ...\n', fileA);
        idx0 = find(fileA == '.');
        idx1 = find(fileA == '/');
        if idx0(end) < idx1(end)
            fprintf('not valid file name; skip!\n');
            continue;
        end
        head = fileA(idx1(end)+1:idx0(end)-1);
        gmw_cmd = '../drivers/gmw';
        % gmw81
        fileL = [ head, '_gmw81_L.mtx' ];
        fileP = [ head, '_gmw81_P.mtx' ];
        fileE = [ head, '_gmw81_E.mtx' ];
        cmd = [ gmw_cmd, ' -gmw81 ', fileA, ' ', fileL, ' -P=', fileP, ' -E=', fileE ];
        fprintf('%s\n', cmd);
        system(cmd);
        % gmw1
        fileL = [ head, '_gmw1_L.mtx' ];
        fileP = [ head, '_gmw1_P.mtx' ];
        fileE = [ head, '_gmw1_E.mtx' ];
        cmd = [ gmw_cmd, ' -gmw1 ', fileA, ' ', fileL, ' -P=', fileP, ' -E=', fileE ];
        fprintf('%s\n', cmd);
        system(cmd);
        % gmw2
        fileL = [ head, '_gmw2_L.mtx' ];
        fileP = [ head, '_gmw2_P.mtx' ];
        fileE = [ head, '_gmw2_E.mtx' ];
        cmd = [ gmw_cmd, ' -gmw2 ', fileA, ' ', fileL, ' -P=', fileP, ' -E=', fileE ];
        fprintf('%s\n', cmd);
        system(cmd);
    end
    aline = fgets(fid);
end

fclose(fid);
