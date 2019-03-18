%Inicialización de Varibles Iniciales
A = [5 -1  1;
     2  8 -1;
    -1  1  4];
b = [10; 11; 3];
I = eye(3);
success = 0;

while(success < 3)
  B = randi([-5, 5], 3, 3); %Generación de matriz B de forma aleatoria
  if (det(B) != 0)          %Confirmar que la matriz sea no singular
    B_inv = inv(B);         %Inversa de B
    G = I - B_inv*A;        %Calculo de matriz G
    eigen = eig(G);         %Calculo de valores eigen
    pG = max(abs(eigen));   %Calculo Radio espectral
    if (pG < 1)             %Verificar que el radio espectral cumpla
       printf("Se encontró la matriz ");
       display(B);          %Mostrar la matriz B encontrada
       printf("El radio espectral es %f \n", pG);
       success = success + 1;
    end
  end
end