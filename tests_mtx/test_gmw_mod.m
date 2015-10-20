fid = fopen('mtx_list.txt');

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
        A = mmread(fileA);
        head = fileA(idx1(end)+1:idx0(end)-1);
        % gmw81
        fileL = [ head, '_gmw81_L.mtx' ];
        fileP = [ head, '_gmw81_P.mtx' ];
        fileE = [ head, '_gmw81_E.mtx' ];
        L = mmread(fileL);
        P = mmread(fileP);
        E = mmread(fileE);
        R = L*L'-P*(A+E)*P';
        lambda = eig(A);
        minE = -min(lambda);
        if minE > 0
            r2 = max(E(:))/minE;
        else
            r2 = max(E(:));
        end
        r2 = full(r2);
        fprintf('gmw81: r2=%.3f, zeta=%d\n', r2, floor(log10(cond(L*L'))));
        fprintf('res=%f\n', norm(R,'fro')/norm(A+E,'fro'));
        % gmw1
        fileL = [ head, '_gmw1_L.mtx' ];
        fileP = [ head, '_gmw1_P.mtx' ];
        fileE = [ head, '_gmw1_E.mtx' ];
        L = mmread(fileL);
        P = mmread(fileP);
        E = mmread(fileE);
        R = L*L'-P*(A+E)*P';
        lambda = eig(A);
        minE = -min(lambda);
        if minE > 0
            r2 = max(E(:))/minE;
        else
            r2 = max(E(:));
        end
        r2 = full(r2);
        fprintf('gmw1: r2=%.3f, zeta=%d\n', r2, floor(log10(cond(L*L'))));
        fprintf('res=%f\n', norm(R,'fro')/norm(A+E,'fro'));
        % gmw2
        fileL = [ head, '_gmw2_L.mtx' ];
        fileP = [ head, '_gmw2_P.mtx' ];
        fileE = [ head, '_gmw2_E.mtx' ];
        L = mmread(fileL);
        P = mmread(fileP);
        E = mmread(fileE);
        R = L*L'-P*(A+E)*P';
        lambda = eig(A);
        minE = -min(lambda);
        if minE > 0
            r2 = max(E(:))/minE;
        else
            r2 = max(E(:));
        end
        r2 = full(r2);
        fprintf('gmw2: r2=%.3f, zeta=%d\n', r2, floor(log10(cond(L*L'))));
        fprintf('res=%f\n', norm(R,'fro')/norm(A+E,'fro'));
    end
    aline = fgets(fid);
end

fclose(fid);
