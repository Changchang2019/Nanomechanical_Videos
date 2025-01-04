function [Lref,davg_mode1,davg_mode2,davg_mode3] = EVM_joel_v6(I,Q)

C = [I(I>0 & Q>0) ; Q(I>0 & Q>0)];

I_select = C(1,:);
Q_select = C(2,:);

% Calculate Lref
I_select_ref = mean(I_select);
Q_select_ref = mean(Q_select);
Lref = (I_select_ref^2 + Q_select_ref^2)^0.5;

% Calculate davg
davg_mode1 = mean(((I_select - I_select_ref).^2 + (Q_select - Q_select_ref).^2).^0.5);
davg_mode2 = (var(I_select(I_select<2*I_select_ref)) + var(Q_select(Q_select<2*Q_select_ref)))^0.5;
davg_mode3 = (mean((I_select - I_select_ref).^2) + mean((Q_select - Q_select_ref).^2)).^0.5;
end

