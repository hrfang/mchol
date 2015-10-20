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
        se_cmd = '../drivers/se';
        % se90
        fileL = [ head, '_se90_L.mtx' ];
        fileP = [ head, '_se90_P.mtx' ];
        fileE = [ head, '_se90_E.mtx' ];
        cmd = [ se_cmd, ' -se90 ', fileA, ' ', fileL, ' -P=', fileP, ' -E=', fileE ];
        fprintf('%s\n', cmd);
        system(cmd);
        % se99
        fileL = [ head, '_se99_L.mtx' ];
        fileP = [ head, '_se99_P.mtx' ];
        fileE = [ head, '_se99_E.mtx' ];
        cmd = [ se_cmd, ' -se99 ', fileA, ' ', fileL, ' -P=', fileP, ' -E=', fileE ];
        fprintf('%s\n', cmd);
        system(cmd);
        % se1
        fileL = [ head, '_se1_L.mtx' ];
        fileP = [ head, '_se1_P.mtx' ];
        fileE = [ head, '_se1_E.mtx' ];
        cmd = [ se_cmd, ' -se1 ', fileA, ' ', fileL, ' -P=', fileP, ' -E=', fileE ];
        fprintf('%s\n', cmd);
        system(cmd);
    end
    aline = fgets(fid);
end

fclose(fid);
