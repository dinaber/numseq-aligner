function [alM, template] = my_aligner (seqs)
% Input : matrix with numbers. Each row is treated as seperate number sequence. 
% Output: aligned matrix with numbers. Each coloumn represents a single
% number. zeros stand for space holders.

% Method: 
% The matrix is transffered into adjacency direct matrix. 
% A topological ordering algorithm for directed graphs is
% applied to produce the ordered template including all numbers. 
% The sequences are then mapped to the template.

% Important: 
% ###The zeros are disregarded. 
% ###Any duplication are treated but not triplications!
% ###The algorithem must received acyclic graph - > self loops must be removed.
seqs = abs(seqs);
max_sp = max(seqs(:));
dict = []; count = 1;
covM = zeros(size(seqs,1),size(seqs,2)); %conversion matrix
toremove = [];

%------%------%------%------%------%------%-----%-----%
%%% First checking the sequences and treatment of duplications:

for i = 1:size(seqs,1)
    seq = seqs(i,:);  seq  = seq(seq>0);
    usq = unique(seq);
    if length(usq) < length(seq) %There are duplications
        hist = histc(seq,usq);
        if any(hist>2) %multiplications
            toremove(end+1) = i;
        end
        dup  = usq(hist>1);
        for d = 1:length(dup)
            n = find(seq==dup(d),2);
            if any(any(dict == dup(d))) %already been replaced
                covM(i,n(end)) = dict(dict(:,1)==dup(d),2) - dup(d);
            else
                dict(end+1,1)  = dup(d);
                dict(end,2)    = max_sp+count;
                covM(i,n(end)) = max_sp+count - dup(d);
                count = count+1;
            end
        end
    end
end
cseqs = seqs +covM ;

if ~isempty(toremove)
    sprintf('Remove row %d due to multiplication\n', toremove)
    cseqs(toremove,:) = [];
end
%------%------%------%------%------%------%-----%-----%


[a1,v1,~] = my_adjacency(cseqs);

if graphisdag(sparse(a1))               %graph is acyclic and valid
    order = graphtopoorder(sparse(a1)); %topological ordering of graph
else
    error('Graph is cyclic / clashing order')
end
of = findall_orderFlip(cseqs);

names        = cellfun(@(x) str2double(x), v1);
template     = names(order);

alM = zeros(size(cseqs,1), length(template));

for s = 1:size(cseqs,1)
    cseq     = cseqs(s,:); cseq = cseq(cseq>0); %seq after treatment
    seq      = seqs(s,:);  seq  = seq(seq>0);   %original seq
    temM     = repmat(template, length(cseq),1);
    [~, inx] = find(temM==cseq');
    
    alM(s,inx) = seq;
end
end
