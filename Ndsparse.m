classdef Ndsparse < handle
    
    properties
        entries
        d
        shape
    end
    
    methods
        
        function this = Ndsparse(varargin)
            % Constructor
            
            this.entries = containers.Map('KeyType','uint64','ValueType','double');
            
            % Blank Ndsparse
            if nargin == 0
                this.d = 0;
                this.shape = [];
            
            % Copy other Ndsparse
            elseif isa(varargin{1},'Ndsparse')
                origNdsparse = varargin{1};
                this.entries = copyMap(origNdsparse.entries);
                this.d = origNdsparse.d;
                this.shape = origNdsparse.shape;
            
            % Single scalar
            elseif length(varargin{1}) == 1
                scalarEntry = varargin{1};
                this.entries(0) = scalarEntry;
                this.d = 0;
                this.shape = [];
            
            % a x N+1 array, where cols 1...N are indices for N-dim matrix
            % and N+1 col is the val and there are a entries
            elseif nargin > 1 && isa(varargin{1},'double')
                listEntries = varargin{1};
                for i = 1:size(listEntries,1)
                    this.addEntry(listEntries(i,1:end-1),listEntries(i,end));
                end
                this.d = size(listEntries,2)-1;
                this.shape = varargin{2};
            
            % N-dim matrix
            else
                denseEntries = varargin{1};
                this.shape = size(denseEntries);
                this.d = length(this.shape);
                for i = 1:prod(this.shape)
                    pos = ind2subs(this.shape,i);
                    val = denseEntries(i);
                    this.addEntry(pos,val);
                end
                
            end
            
            this.removeZeros();
            
        end
        
        function disp(this)
            rep = [num2str(this.d) '-d sparse tensor with ' num2str(this.nnz()) ' nonzero entries\n'];
            indList = sort(cell2mat(keys(this.entries)));
            for i = indList
                rep = [rep '(' num2str(ind2subs(this.shape,i)) ')\t' num2str(this.entries(i)) '\n'];
            end
            fprintf(rep)
        end
        
        function addEntry(this, pos, val)
            % Add entry, adding to element if already present
            ind = subs2ind(this.shape,pos);
            if isKey(this.entries,ind)
                this.entries(ind) = this.entries(ind) + val;
            else
                this.entries(ind) = val;
            end
        end
        
        function removeEntry(this, pos)
            ind = subs2ind(this.shape,pos);
            remove(this.entries,ind);
        end
        
        function numEntries = nnz(this)
            numEntries = this.entries.length;
        end
        
        function removeZeros(this)
            indList = cell2mat(keys(this.entries));
            for i = 1:this.entries.length
                if this.entries(indList(i)) == 0
                    this.removeEntry(ind2subs(this.shape,i));
                end
            end
        end
        
        function prod = ttt(this, spec1, that, spec2)
            % Compute the tensor contraction of this with that, contracting
            %   (inner product) along the corresponding negative indices in
            %   spec1 and spec2 and arranging the product according to the
            %   positive indices in spec1 and spec2. Like Einstein notation.
            
            % Validate dimensions
            assert(this.d == length(spec1), 'Ndsparse:ttt:Tensor1DimMismatch', 'Ndims of 1st tensor (%d) doesn''t match specified dims (%d).', this.d, length(spec1))
            assert(that.d == length(spec2), 'Ndsparse:ttt:Tensor2DimMismatch', 'Ndims of 2nd tensor (%d) doesn''t match specified dims (%d).', that.d, length(spec2))
            
            % Validate inner product dims
            inDims1 = spec1(spec1 < 0);
            inDims2 = spec2(spec2 < 0);
            if ~all(ismember(inDims1, inDims2))
                conDiff1 = setdiff(inDims1, inDims2); % in case there are multiple extra dims
                error('Ndsparse:ttt:Tensor1ExtraInnerProd', '1st tensor has extra inner prod dim (%d)', conDiff1(1))
            end
            if ~all(ismember(inDims2, inDims1))
                conDiff2 = setdiff(inDims2, inDims1);
                error('Ndsparse:ttt:Tensor2ExtraInnerProd', '2nd tensor has extra inner prod dim (%d)', conDiff2(1))
            end
            
            % Validate outer product dims
            out1 = spec1 > 0;
            out2 = spec2 > 0;
            outDims1 = spec1(out1);
            outDims2 = spec2(out2);
            outDims = [outDims1, outDims2];
            assert(length(outDims) == length(unique(outDims)), 'Ndsparse:ttt:RepeatedOuterProdDims', 'Outer product dimensions must be unique')
            assert(all(sort(outDims) == 1:length(outDims)), 'Ndsparse:ttt:InvalidOuterProdDimOrder', 'Outer product dimensions must contain integers 1:nLeftoverDims')
            
            % Make mapping between this's inner product dims with that's inner product dims
            nInDims = length(inDims1);
            [~, in1] = sort(spec1);
            in1 = in1(1:nInDims);
            [~, in2] = sort(spec2);
            in2 = in2(1:nInDims);
            
            % Validate that corresponding inner product shapes match
            for i = 1:nInDims
                assert(this.shape(in1(i)) == that.shape(in2(i)), 'Ndsparse:ttt:InnerProdShapeMismatch', '%d-th inner product shape mismatch. 1st tensor has %d elements and 2nd tensor has %d elements %d', i, this.shape(in1(i)), that.shape(in2(i)))
            end
            
            % Make ordering of outer product dims and product shape
            out = [out1 out2];
            allShape = [this.shape that.shape];
            outShape = allShape(out);
            prodShape = outShape(outDims);
            
            % Initialize product tensor
            prod = Ndsparse();
            prod.shape = prodShape;
            prod.d = length(prod.shape);
            
            % Collect linear indices
            indList1 = cell2mat(keys(this.entries));
            indList2 = cell2mat(keys(that.entries));
            
            % Perform contraction, only taking nonzero entries
            for i = indList1
                
                pos1 = ind2subs(this.shape, i);
                inInds1 = pos1(in1);
                outInds1 = pos1(out1);
                
                for j = indList2
                    
                    pos2 = ind2subs(that.shape,j);
                    inInds2 = pos2(in2);
                    outInds2 = pos2(out2);
                    
                    if all(inInds1 == inInds2)
                        outInds = [outInds1 outInds2];
                        pos = outInds(outDims);
                        val = this.entries(i) * that.entries(j);
                        prod.addEntry(pos, val);
                    end
                    
                end
            end
            
        end
        
        function transpose(this, permutation)
            % Transpose/permute dimensions according to permutation.
            % Note: there's probably a cleaner way to do this - inplace
            
            newMap = containers.Map('KeyType','uint64','ValueType','double');
            newShape = this.shape(permutation);
            
            indList = cell2mat(keys(this.entries));
            for i = 1:this.entries.length
                subsOld = ind2subs(this.shape,indList(i));
                subsNew = subsOld(permutation);
                indNew = subs2ind(newShape,subsNew);
                newMap(indNew) = this.entries(indList(i));
            end
            
            this.entries = newMap;
            this.shape = newShape;
        end
        
        function out = full(this)
            % Convert Ndsparse tensor to Matlab full N-dim matrix
            out = zeros(this.shape);
            inds = cell2mat(keys(this.entries));
            for i = inds
                out(i) = this.entries(i); % linear indexing
            end
        end
        
        
    end
    
end