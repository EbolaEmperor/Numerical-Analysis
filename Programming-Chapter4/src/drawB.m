A = [0 -3.5; 0 -3; 0 -2.5; 0 -2; 0 -1.75; 0 -1.5; 0 -1.25; 0 -1; 0 -0.875; 0 -0.75; 0 -0.625; 0 -0.5; 0 0; 
     0 0.5; 0 0.625; 0 0.75; 0 0.875; 0 1; 0 1.25; 0 1.5; 0 1.75; 0 2; 0 2.5; 0 3; 0 3.5];

axis([-4 4 -0.1 0.1]);
hold on
scatter(A(:,2),A(:,1),'filled');
figure;

B = [0 -0.375; 0 -0.25; 0 -0.125; 0 0; 0 0.125; 0 0.25; 0 0.375];

axis([-4 4 -0.1 0.1]);
hold on
scatter(A(:,2),A(:,1),'filled');
scatter(B(:,2),B(:,1),'filled');
