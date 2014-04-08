function rc = FMGrestrict2D(rf)
% removes half of the interior grid points

temp1 = 1:length(rf);
temp2 = temp1';
temp1(:,1:2:end)=[];
temp2(1:2:end,:)=[];
rc = rf(temp1,:);
rc = rc(:,temp2);

end