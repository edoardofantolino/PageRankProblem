close all
clear
clc

prompt = 'Which requirement do you want to see? ';
req = input(prompt,'s');

data = dlmread('edges_file_matlab.txt');
G = sparse(data(:, 1), data(:, 2), 1);
G = [G zeros(length(G),1)];

r_j = sum(G,2);      % indegree
c_j = sum(G,1);       % outdegree

switch req
    case '1'
        % (*1) G SPARSITY PATTERN
        figure(1)
        spy(G)
        title('G matrix sparsity pattern')
        
    case '2'
        % (*2)
        disp('See the report for further details')
        
    case '3'
        % (*3) INDEGREEE OUTDEGREE
        r_j      % indegree
        c_j      % outdegree
        
    case '4'
        % (**4)
        disp('See the report for further details')
        
    case '5'
        % (*5) HYPERLING MATRIX M
        
        M = G(:,1)/c_j(1);
        for k=2:length(G)
            if c_j(k)~=0 
                M = [M G(:,k)/c_j(k)];
            else 
                M = [M zeros(length(G), 1)];
            end
        end
        
        figure(1)
        spy(G)
        title('G matrix sparsity pattern')
        
        figure(2)
        spy(M)
        title('M matrix sparsity pattern')
        
    case '6'
        % (*6)
        disp('See the report for further details')
        
    case '7'
        % (*7)
        Mt = G(:,1)/c_j(1);
        for k=2:length(G)
            if c_j(k)~=0 
                Mt = [Mt G(:,k)/c_j(k)];
            else 
                Mt = [Mt (1/length(G))*ones(length(G),1)];
            end
        end
        
        figure(3)
        spy(Mt)
        title('Mt matrix sparsity pattern')
        
    case '8'
        % (*8)
        disp('See the report for further details')
        
    case '9'
        % (**9)      
        disp('See the report for further details')
        
    case '10'
        % (*10)
        Mt = G(:,1)/c_j(1);
        for k=2:length(G)
            if c_j(k)~=0 
                Mt = [Mt G(:,k)/c_j(k)];
            else 
                Mt = [Mt (1/length(G))*ones(length(G),1)];
            end
        end
        
        lMATLAB = eigs(Mt);
        [lambda, x, iter] = power_method(Mt);       
        rel_err = abs(lMATLAB(1) - lambda)/ abs(lMATLAB(1));
        disp(['labmda power method: ', num2str(lambda)])
        disp(['labmda eigs: ', num2str(lMATLAB(1))])
        disp(['The relative error between lambda calculated with ', ...
            'eigs() and the power_method() is: ', num2str(rel_err)])
        
    case '11'
        % (*11)
        alpha = 0.85;
        delta = (1-alpha)/length(G);
        A = alpha*G(:,1)/c_j(1)+delta;
        for k=2:length(G)
            if c_j(k)~=0 
                A = [A alpha*G(:,k)/c_j(k)+delta];
            else
                A = [A (1/length(G))*ones(length(G),1)];
            end
        end
        
        sum(A)
        disp(['number of elements: ', num2str(length(A)*length(A))])
        disp(['number of nonzeros: ', num2str(nnz(A))])
        
        
    case '12'
        % (*12)
        alpha = 0.85;
        delta = (1-alpha)/length(G);
        A = alpha*G(:,1)/c_j(1)+delta;
        for k=2:length(G)
            if c_j(k)~=0 
                A = [A alpha*G(:,k)/c_j(k)+delta];
            else
                A = [A (1/length(G))*ones(length(G),1)];
            end
        end
        
    case '13'
        % (**13)
        disp('See the report for further details')
        
    case '14'
        % (**14)
        alpha = 0.85;
        delta = (1-alpha)/length(G);
        A = alpha*G(:,1)/c_j(1)+delta;
        for k=2:length(G)
            if c_j(k)~=0 
                A = [A alpha*G(:,k)/c_j(k)+delta];
            else
                A = [A (1/length(G))*ones(length(G),1)];
            end
        end
        
        [V,L] = eigs(A,1);
        xMATLAB = V/sum(V);
        [lambda, x, iter] = power_method(A);
        x = x/sum(x);
        
        lambda_rel_err = abs(L-lambda) / abs(L);
        disp(['The relative error between lambda calculated with ', ...
            'eigs() and the power_method() is: ', num2str(lambda_rel_err)])
        
        x_rel_err = norm(xMATLAB - x) / norm(xMATLAB);
        disp(['The relative error between eigenvector calculated with ', ...
            'eigs() and the power_method() is: ', num2str(x_rel_err)])
        
    case '15'
        % (**15) and (***16)
        alpha = 0.85;
        delta = (1-alpha)/length(G);
        e = ones(length(G),1);
        z = (1/length(G))*ones(length(G),1);

        d_j = zeros(length(G),1);
        for k=1:length(G)
            if c_j(k) ~= 0
                d_j(k) = 1/c_j(k);
                z(k) = delta;
            end
        end

        D = spdiags(d_j,0,length(G),length(G));
        A = alpha*G*D + e*z';
        
        [X,L]=eigs(A,1);
%         disp('PageRank matlab: ')
        X = X'/sum(X,1);
        
        [lambda,x,iter] = power_method(A);
        disp('PageRank power method: ')
        x = x'/sum(x,1);
        x_rel_err = norm(X - x) / norm(X);
        disp(['The relative error between eigenvector calculated with ', ...
            'eigs() and the power_method() is: ', num2str(x_rel_err)])

        [lambda2,x2,iter2] = sparse_power_method(alpha,G,D,e,z);
        disp('PageRank new sparse power method: ')
        x2 = x2'/sum(x2,1);
        x_rel_err = norm(X - x2) / norm(X);
        disp(['The relative error between eigenvector calculated with ', ...
            'eigs() and the sparse_power_method() is: ', num2str(x_rel_err)])
        
        [lambda3,x3,iter3] = no_mm_power_method(alpha,G,D,e,z, c_j);
        disp('PageRank no multiplication matrix vector power method: ')
        x3 = x3'/sum(x3,1);
        x_rel_err = norm(X - x3) / norm(X);
        disp(['The relative error between eigenvector calculated with ', ...
            'eigs() and the no_mm_sparse_power_method() is: ', num2str(x_rel_err)])
        
    case '16'
        % (***16)
        disp('Case 16 is in case 15, no_mm_power_method')
        
    case '17'
        % (**17)
        Mt = G(:,1)/c_j(1);
        for k=2:length(G)
            if c_j(k)~=0 
                Mt = [Mt G(:,k)/c_j(k)];
            else 
                Mt = [Mt (1/length(G))*ones(length(G),1)];
            end
        end
        
        [V,L] = eigs(Mt,2);
        [lambda, x, iter] = power_method(Mt);
        [lambda_2, x_2] = power_deflation(Mt,lambda,x);
        
        rel_err = abs(L(2,2) - lambda_2) / abs(L(2,2));
        disp(['The absolute value of lambda_2 should be almost 1 ---> ', num2str(lambda_2)])
        disp(['The relative error between lambda2 calculated with ', ...
            'eigs() and the power_method() is: ', num2str(rel_err)])
        
    case '18'
        % (**18)
        disp('See the report for further details')        
    case '19'
        % (*19)
        tic;
        alpha = 0.85;
        delta = (1-alpha)/length(G);
        e = ones(length(G),1);
        z = (1/length(G))*ones(length(G),1);

        d_j = zeros(length(G),1);
        for k=1:length(G)
            if c_j(k) ~= 0
                d_j(k) = 1/c_j(k);
                z(k) = delta;
            end
        end

        D = spdiags(d_j,0,length(G),length(G));
        [lambda, x, iter] = power_method(alpha*G*D+e*z');
        toc;
        
        x = x/sum(x);
        N = alpha*G*D+e*z';
                
        bar(x, 'BarWidth', 4)
        title('PageRank bar graph')
        
        dozen = 12;
        
        [value, page] = maxk(x, dozen);
        
        for k=1:dozen
            disp(['rank: ', num2str(k), ...
                  ', page number: ', num2str(page(k)), ...
                  ', value rank: ', num2str(value(k)), ...
                  ', outdegree: ', num2str(c_j(page(k))), ...
                  ', indegree: ', num2str(r_j(page(k)))])
        end
        
        all_features = ["187", "1", "188", "1030", "428", ...
                        "76", "101", "382", "186", "1320", ...
                        "1319", "385"];
        figure(2)
        bar(value, 'k') 
        grid on
        set(gca,'xticklabel',all_features)

        
%         standings = maxk(x, dozen);
%         pages = zeros(dozen, 2);
% 
%         for k=1:dozen
%             k;
%             values = find(x==standings(k));
%             if length(values) == 1
%                 pages(k,1) = values;
%                 pages(k,2) = 0;
%             else
%                 pages(k,:) = values;
%             end
%         end
% 
%         for k=1:dozen
%             disp(['page number: ', num2str(pages(k)), ...
%                   ', outdegree: ', num2str(c_j(pages(k))), ...
%                   ', indegree: ', num2str(r_j(pages(k)))])
%         end
    case '20'
        disp('We are in 20')
        tic;
        alpha = 0.85;
        delta = (1-alpha)/length(G);
        e = ones(length(G),1);
        z = (1/length(G))*ones(length(G),1);

        d_j = zeros(length(G),1);
        for k=1:length(G)
            if c_j(k) ~= 0
                d_j(k) = 1/c_j(k);
                z(k) = delta;
            end
        end

        D = spdiags(d_j,0,length(G),length(G));
        [lambda, x, iter] = power_method(alpha*G*D+e*z');
        toc;
        
        
        x = x/sum(x);
        A = alpha*G*D+e*z';
        
        [V,L] = eigs(A);
        vmat = V(:,1);
        
        v = rand(length(A),1);
        v = v/norm(v);
        vstart = v;

        num_iter = 100;
        errors = zeros(num_iter,1);
        vseq = zeros(length(A), num_iter);
        for i=1:num_iter
           v = A*v; 
           v = v/norm(v);
           vseq(:,i) = v;
           errors(i) = norm(vmat-v)/norm(vmat);
           i
        end
        
        figure(1)
        semilogy(errors, 'o')
        grid on
        
        figure(2)
        plot(vseq, 'ko')
        hold on 
        grid on
        scatter([1:2500], vseq(:,1), 'ro', 'filled')
        scatter([1:2500], vmat, 'go', 'filled')
        
end
