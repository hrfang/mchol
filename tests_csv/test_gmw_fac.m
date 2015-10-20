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
        A = csvread(fileA);
        % gmw81
        [L, P, E] = gmw81(A);
        head = fileA(idx1(end)+1:idx0(end)-1);
        fileL = [ head, '_gmw81_L.csv' ];
        fileP = [ head, '_gmw81_P.csv' ];
        fileE = [ head, '_gmw81_E.csv' ];
        csvwrite(fileL, L);
        csvwrite(fileP, full(P));
        csvwrite(fileE, full(E));
        % gmw1
        [L, P, E] = gmw1(A);
        fileL = [ head, '_gmw1_L.csv' ];
        fileP = [ head, '_gmw1_P.csv' ];
        fileE = [ head, '_gmw1_E.csv' ];
        csvwrite(fileL, L);
        csvwrite(fileP, full(P));
        csvwrite(fileE, full(E));
        % gmw2
        [L, P, E] = gmw2(A);
        fileL = [ head, '_gmw2_L.csv' ];
        fileP = [ head, '_gmw2_P.csv' ];
        fileE = [ head, '_gmw2_E.csv' ];
        csvwrite(fileL, L);
        csvwrite(fileP, full(P));
        csvwrite(fileE, full(E));
    end
    aline = fgets(fid);
end

fclose(fid);
