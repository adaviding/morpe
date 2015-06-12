ddp = 0.01;
dPrime = 0:ddp:5;
h = zeros(size(dPrime));
for i=1:length(dPrime)
	h(i) = Mcl_Hd(dPrime(i));
end

figure;
plot(dPrime, h, 'Color', [0 0 0], 'LineWidth', 6);
hold on;
grid on;
%dp = dPrime/1.1595;
%plot(dPrime, (1./(1+dp)), 'Color', [0 0 1], 'LineWidth', 2);
%plot(dPrime, (1./(1+dp.^2)), 'Color', [0 0.75 0], 'LineWidth', 2);
%plot(dPrime, (1./(1+dp.^max(1,dp.^0.78))), 'Color', [1 0 0], 'LineWidth', 4);
%plot(dPrime, exp(-dPrime*log(2)), 'Color', [1 0 0], 'LineWidth', 3);
xlabel('dPrime');
ylabel('Conditional Entropy');
