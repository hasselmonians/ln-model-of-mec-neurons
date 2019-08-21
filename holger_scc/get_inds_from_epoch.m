function [indices] = get_inds_from_epoch(root,epochs)

root.epoch = [-inf inf];
inds = cell(1,size(epochs,1));
for i = 1:size(epochs,1)
    inds{1,i} = [find(abs(root.ts-epochs(i,1)) == min(abs(root.ts-epochs(i,1)))):find(abs(root.ts-epochs(i,2)) == min(abs(root.ts-epochs(i,2))))];
end

indices = cell2mat(inds);
