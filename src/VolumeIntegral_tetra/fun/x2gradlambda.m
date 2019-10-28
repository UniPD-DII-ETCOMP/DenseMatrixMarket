function [gradlambda] = x2gradlambda(norm,area,vol)
gradlambda(1:3,1:4)=0;
vol3=3.0d0*vol;
for i=1:4
gradlambda(:,i)=-norm(:,i)*area(i)/vol3;
end
end 