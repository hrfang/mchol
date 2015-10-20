path(path, '../mex/');
fid = fopen('csv_list.txt');

aline = fgetl(fid);
while ischar(aline)
    fileA = strtrim(aline);  % trim leading and tailing spaces etc.
    if length(fileA)
        fprintf('# processing %s ...\n', fileA);
        idx0 = find(fileA == '.');
        idx1 = find(fileA == '/');
        if idx0(end) < idx1(end)
            fprintf('not valid file name; skip!\n');
            continue;
        end
        A = csvread(fileA);
        head = fileA(idx1(end)+1:idx0(end)-1);
        % se90
        fileL = [ head, '_se90_L.csv' ];
        fileP = [ head, '_se90_P.csv' ];
        fileE = [ head, '_se90_E.csv' ];
        L = csvread(fileL);
        P = csvread(fileP);
        E = csvread(fileE);
        R = L*L'-P*(A+E)*P';
        lambda = eig(A);
        minE = -min(lambda);
        if minE > 0
            r2 = max(E(:))/minE;
        else
            r2 = max(E(:));
        end
        r2 = full(r2);
        fprintf('se90: r2=%.3f, zeta=%d\n', r2, floor(log10(cond(L*L'))));
        fprintf('res=%f\n', norm(R,'fro')/norm(A+E,'fro'));
        % se99
        fileL = [ head, '_se99_L.csv' ];
        fileP = [ head, '_se99_P.csv' ];
        fileE = [ head, '_se99_E.csv' ];
        L = csvread(fileL);
        P = csvread(fileP);
        E = csvread(fileE);
        R = L*L'-P*(A+E)*P';
        lambda = eig(A);
        minE = -min(lambda);
        if minE > 0
            r2 = max(E(:))/minE;
        else
            r2 = max(E(:));
        end
        r2 = full(r2);
        fprintf('se99: r2=%.3f, zeta=%d\n', r2, floor(log10(cond(L*L'))));
        fprintf('res=%f\n', norm(R,'fro')/norm(A+E,'fro'));
        % se1
        fileL = [ head, '_se1_L.csv' ];
        fileP = [ head, '_se1_P.csv' ];
        fileE = [ head, '_se1_E.csv' ];
        L = csvread(fileL);
        P = csvread(fileP);
        E = csvread(fileE);
        R = L*L'-P*(A+E)*P';
        lambda = eig(A);
        minE = -min(lambda);
        if minE > 0
            r2 = max(E(:))/minE;
        else
            r2 = max(E(:));
        end
        r2 = full(r2);
        fprintf('se1: r2=%.3f, zeta=%d\n', r2, floor(log10(cond(L*L'))));
        fprintf('res=%f\n', norm(R,'fro')/norm(A+E,'fro'));
    end
    aline = fgets(fid);
end

fclose(fid);
