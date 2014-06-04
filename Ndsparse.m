classdef Ndsparse < handle
    
    properties
        entries
        d
        shape
    end
    
    methods
        
        function self = Ndsparse(varargin)
            % Constructor
            
            self.entries = containers.Map('KeyType','uint64','ValueType','double');
            
            % Blank Ndsparse
            if nargin == 0
                self.d = 0;
                self.shape = [];
            
            % Copy other Ndsparse
            elseif isa(varargin{1},'Ndsparse')
                origNdsparse = varargin{1};
                self.entries = copyMap(origNdsparse.entries);
                self.d = origNdsparse.d;
                self.shape = origNdsparse.shape;
            
            % Single scalar
            elseif length(varargin{1}) == 1
                scalarEntry = varargin{1};
                self.entries(0) = scalarEntry;
                self.d = 0;
                self.shape = [];
            
            % a x N+1 array, where cols 1...N are indices for N-dim matrix
            % and N+1 col is the val and there are a entries
            elseif nargin > 1 && isa(varargin{1},'double')
                listEntries = varargin{1};
                for i = 1:size(listEntries,1)
                    self.addEntry(listEntries(i,1:end-1),listEntries(i,end));
                end
                self.d = size(listEntries,2)-1;
                self.shape = varargin{2};
            
            % N-dim matrix
            else
                denseEntries = varargin{1};
                self.shape = size(denseEntries);
                self.d = length(self.shape);
                for i = 1:prod(self.shape)
                    pos = ind2subs(self.shape,i);
                    val = denseEntries(i);
                    self.addEntry(pos,val);
                end
                
            end
            
            self.removeZeros();
            
        end
        
        function disp(self)
            rep = [num2str(self.d) '-d sparse tensor with ' num2str(self.nnz()) ' nonzero entries\n'];
            indList = sort(cell2mat(keys(self.entries)));
            for i = indList
                rep = [rep '(' num2str(ind2subs(self.shape,i)) ')\t' num2str(self.entries(i)) '\n'];
            end
            fprintf(rep)
        end
        
        function addEntry(self,pos,val)
            % Add entry, adding to element if already present
            ind = subs2ind(self.shape,pos);
            if isKey(self.entries,ind)
                self.entries(ind) = self.entries(ind) + val;
            else
                self.entries(ind) = val;
            end
        end
        
        function removeEntry(self,pos)
            ind = subs2ind(self.shape,pos);
            remove(self.entries,ind);
        end
        
        function numEntries = nnz(self)
            numEntries = self.entries.length;
        end
        
        function removeZeros(self)
            indList = cell2mat(keys(self.entries));
            for i = 1:self.entries.length
                if self.entries(indList(i)) == 0
                    self.removeEntry(ind2subs(self.shape,i));
                end
            end
        end
        
        function out = ttt(self,other,varargin)
            
            if nargin < 3
                con1 = [];
                con2 = [];
            elseif nargin < 5
                con1 = varargin{1};
                con2 = varargin{2};
            else
                error('Error: Ndsparse ttt too many arguments')
            end
            keep1 = setdiff(1:self.d,con1);
            keep2 = setdiff(1:other.d,con2);
            
            out = Ndsparse();
            out.shape = [self.shape(keep1) other.shape(keep2)];
            out.d = length(out.shape);
            
            indList1 = cell2mat(keys(self.entries));
            indList2 = cell2mat(keys(other.entries));
            
            for i = indList1
                
                pos1 = ind2subs(self.shape,i);
                conIdx1s = pos1(con1);
                keepIdx1s = pos1(keep1);
                
                for j = indList2
                    
                    pos2 = ind2subs(other.shape,j);
                    conIdx2s = pos2(con2);
                    keepIdx2s = pos2(keep2);
                    
                    if all(conIdx1s == conIdx2s)
                        pos = [keepIdx1s keepIdx2s];
                        val = self.entries(i)*other.entries(j);
                        out.addEntry(pos,val);
                    end
                    
                end
            end
            
        end
        
        function transpose(self,permutation)
            % Note: there's probably a cleaner way to do this
            
            newMap = containers.Map('KeyType','uint64','ValueType','double');
            newShape = self.shape(permutation);
            
            indList = cell2mat(keys(self.entries));
            for i = 1:self.entries.length
                subsOld = ind2subs(self.shape,indList(i));
                subsNew = subsOld(permutation);
                indNew = subs2ind(newShape,subsNew);
                newMap(indNew) = self.entries(indList(i));
            end
            
            self.entries = newMap;
            self.shape = newShape;
        end
        
    end
    
end