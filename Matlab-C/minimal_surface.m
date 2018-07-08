% mex -largeArrayDims mcc.cpp

[x, t] = readObj('Cat_head.obj');

initx=x;
%% draw 2 copies of the image
figure; set(gcf, 'Units', 'normalized', 'Position', [0.05,0.05,.8,.8]);
subplot(121); trimesh(t, x(:,1), x(:,2), x(:,3), 'edgecolor', 'k'); axis off; axis equal; title('input');
subplot(122); h = trimesh(t, x(:,1), x(:,2), x(:,3), 'edgecolor', 'k'); axis off; axis equal; title('output');

%% TODO: find interior vertices
nv = size(x, 1);
MVtx2Vtx = sparse(t, t(:, [2 3 1]), true, nv, nv);% Adjective Matrix£¿£¿

[e1, e2] = find( xor( MVtx2Vtx, MVtx2Vtx' ) );
B = unique([e1; e2]);  % index of boundary vertices//  Q
I = setdiff(1:nv, B);% NOT-boundary


runLocalVersion = true;

if runLocalVersion
   
    %% compute minimal surface using local approach
%    for it=1:10
        L = laplacian(x, t);
        
        %L = L./full( diag(L) );%%%?
        L = L./diag(L);
        
        L = spdiags(zeros(nv,1), 0, L);%%%   Q
    
        L(B,:)=zeros(size(B,1),nv);
        
        for itr=1:size(B,1)
            L(B(itr),B(itr))=-1;   
        end

        FL=full(L);
        
% for it=1:10
%         while true
%             xI = -L*x;
%              if norm(x - xI)<1e-5
%                   fprintf('converged');
%                  break;
%              end
% 
%             x=xI;
%             y=x;
% 
%             set(h, 'Vertices', x); %drawnow; pause(0.1);%%display the process
%         end
% end
tic
xI=mcc(L,x);
toc
set(h, 'Vertices', xI); 
else
    %% compute minimal surface using global approach
    for it=1:10
        L = laplacian(x, t);
        x(I,:) = -L(I,I)\(L(I,B)*x(B,:));%%
        set(h, 'Vertices', x); drawnow; pause(0.1);
    end
end
