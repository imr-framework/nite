function Tp = label2temp(label, Tp_label_vec)

Tp = zeros(size(label));
for npt = 1:length(Tp)
    Tp(npt) = Tp_label_vec(label(npt)+1); %label goes from 0 to N-1
end
