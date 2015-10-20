path(path, '../mex/');
fid = fopen('csv_list.txt');

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
        A = csvread(fileA);
        % se90
        [L, P, E] = se90(A);
        fileL = [ head, '_se90_L.csv' ];
        fileP = [ head, '_se90_P.csv' ];
        fileE = [ head, '_se90_E.csv' ];
        csvwrite(fileL, L);
        csvwrite(fileP, full(P));
        csvwrite(fileE, full(E));
        % se99
        [L, P, E] = se99(A);
        fileL = [ head, '_se99_L.csv' ];
        fileP = [ head, '_se99_P.csv' ];
        fileE = [ head, '_se99_E.csv' ];
        csvwrite(fileL, L);
        csvwrite(fileP, full(P));
        csvwrite(fileE, full(E));
        % se1
        [L, P, E] = se1(A);
        fileL = [ head, '_se1_L.csv' ];
        fileP = [ head, '_se1_P.csv' ];
        fileE = [ head, '_se1_E.csv' ];
        csvwrite(fileL, L);
        csvwrite(fileP, full(P));
        csvwrite(fileE, full(E));
    end
    aline = fgets(fid);
end

fclose(fid);
