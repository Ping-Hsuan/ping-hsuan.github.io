**December 06, 2023** - Note 1

## Generate the POD basis functions using SVD approach

Given a snapshot matrix $X \in \mathbb{R}^{\mathcal{N} \times K}$, where $\mathcal{N}$ is the degrees-of-freedom and $K$ is the number of snapshots:

1. Proper orthogonal decomposition (POD)
2. Singular value decomposition (SVD)

### $l_2$ Case
#### POD
1. $G_{l_2} = X^T X$
2. $G_{l_2} = Q \Lambda Q^{-1}$, where $Q$ is the eigenvector matrix and $\Lambda$ is the eigenvalue matrix
3. Reduced basis matrix $\Xi := [ \xi_1 \cdots \xi_K] = X Q$ (**The reduced basis function $\xi_i$ is $l^2$-orthogonal**)

```
G=zeros(ndump,ndump);
for j=1:ndump
    disp(['gengram ',num2str(j)]);
    Buj=(us0(:,j));
    for i=1:ndump
        G(i,j)=dot(us0(:,i),Buj);
    end 
end

G=.5*(G+G'); [V,D]=eig(G);
D=flip(diag(D)); V=flip(V,2);

ub=zeros((Nx+1)*(Ny+1)*2,ndump);
ub(:,1)=u0; ub(:,2:end)=us0*V(:,1:ndump-1);

% Normalization
for i=2:ndump
    ub(:,i)=ub(:,i)./sqrt(dot(ub(:,i),Bb(ub(:,i))));
end          
```


#### SVD
1. $X = U\Sigma V^T$
2. $G_{l_2} = X^TX = V\Sigma U^T U \Sigma V^T$ = $V \Sigma^2 V^T$ = $Q \Lambda Q^{-1}$
3. Reduced basis matrix $\Xi = XV = U\Sigma$

```
[U_,S_,V_]=svd(us0,"econ"); % U_ is already normalized
```

### $L^2$ Case
#### POD
1. $G_{L_2} = X^T B X$
2. $G_{L_2} = Q \Lambda Q^{-1}$
3. $\Xi = XQ$  (**The reduced basis function $\xi_i$ is $L^2$-orthogonal**)

```
G=zeros(ndump,ndump);
for j=1:ndump
    disp(['gengram ',num2str(j)]);
    Buj=Bu(us0(:,j));
    for i=1:ndump
        G(i,j)=dot(us0(:,i),Buj);
    end 
end

G=.5*(G+G'); [V,D]=eig(G);
D=flip(diag(D)); V=flip(V,2);

ub=zeros((Nx+1)*(Ny+1)*2,ndump);
ub(:,1)=u0; ub(:,2:end)=us0*V(:,1:ndump-1);

% Normalization
for i=2:ndump
    ub(:,i)=ub(:,i)./sqrt(dot(ub(:,i),Bb(ub(:,i))));
end          
```

#### SVD
1. $B^{1/2} X = U\Sigma V^T$
2. $G_{L_2} = X^TB X = V\Sigma U^T U \Sigma V^T$ = $V \Sigma^2 V^T$ = $Q \Lambda Q^{-1}$
3. $\Xi = X V = B^{-1/2} U \Sigma$
```
[U_,S_,V_]=svd(a2u(Bby.^(1/2),Bbx.^(1/2),us0),"econ");

% Need to take B^{-1/2} to get the POD basis
p1(a2u(inv(Bby.^(1/2)),inv(Bbx.^(1/2)),U_(:,idump))); title(['Velocity Magnitude  ' cccc]);
```

<script type="text/javascript"
      src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML">
    </script>
