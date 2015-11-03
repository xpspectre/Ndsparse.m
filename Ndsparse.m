classdef Ndsparse < handle
    
    properties
        entries
        d
        shape
    end
    
    methods
        
        function this = Ndsparse(varargin)
            % Construct Ndsparse object
            % 
            % Syntax:
            %   this = Ndsparse()
            %   this = Ndsparse(that)
            %   this = Ndsparse(scalarEntry)
            %   this = Ndsparse(listEntries, shape)
            %   this = Ndsparse(denseEntries)
            %
            % Description:
            %   this = Ndsparse() initializes a new, empty Ndsparse object.
            %   this = Ndsparse(that) copies the existing Ndsparse that.
            %   this = Ndsparse(scalarEntry) creates a new Ndsparse with
            %       dimension 0, shape [], and single value scalarEntry at index 0.
            %   this = Ndsparse(listEntries, shape) creates a new Ndsparse using
            %       listEntries. listEntries is an a x N+1 double matrix with cols
            %       1...N as indices in N dims and col N+1 as the value. shape
            %       is a 1 x d vector of nonnegative integers indicating the
            %       sizes of the dimensions, which is needed because the sparse
            %       entries don't otherwise give the shape of the bounding
            %       matrix.
            %   this = Ndsparse(denseEntries) creates a new Ndsparse from the
            %       dense, N-dim matrix denseEntries
            %
            % Ndsparse fields:
            %   entries [ Map ]
            %       The nonzero entries (can contain explicit 0's in some cases)
            %       of the sparse matrix, with key/value pairs:
            %       keys [ char ]
            %           Position indices represented by a string of integers
            %           separated by spaces. The string representation is
            %           obtained from num2char of a 1 x d vector of nonnegative
            %           integers.
            %       values [ double ]
            %           Value at key
            %   d [ nonegative scalar integer ] 
            %       Number of dimensions
            %   shape [ 1 x d vector of nonnegative integers ]
            %       Size of each dimension
            
            this.entries = containers.Map('KeyType', 'char', 'ValueType', 'double');
            
            % Blank Ndsparse
            if nargin == 0
                this.d = 0;
                this.shape = [];
            
            % Copy other Ndsparse
            elseif isa(varargin{1}, 'Ndsparse')
                origNdsparse = varargin{1};
                this.entries = copyMap(origNdsparse.entries);
                this.d       = origNdsparse.d;
                this.shape   = origNdsparse.shape;
            
            % Single scalar
            elseif numel(varargin{1}) == 1 && isnumeric(varargin{1})
                scalarEntry = varargin{1};
                this.entries('0') = scalarEntry;
                this.d          = 0;
                this.shape      = [];
            
            % a x N+1 array, where cols 1...N are indices for N-dim matrix
            %   and N+1 col is the val and there are a entries
            elseif nargin == 2 && isnumeric(varargin{1}) && isnumeric(varargin{2})
                listEntries = varargin{1};
                shape       = varargin{2};
                d = size(listEntries,2)-1;
                
                assert(length(shape) == d, 'Ndsparse:Constructor:ShapeDimMismatch', 'Shape dimensions (%d) doesn''t match number of indices (%d)', length(shape), d)
                
                listIndices = listEntries(:,1:d);
                
                assert(isequal(listIndices, floor(listIndices)), 'Ndsparse:Constructor:NonintegerIndices', 'Indices must be integers')
                assert(all(all(listIndices >= 0)), 'Ndsparse:Constructor:NegativeIndices', 'Indices must be nonegative')
                for i = 1:d
                    assert(all(listIndices(:,i) <= shape(i)), 'Ndsparse:Constructor:TooLargeIndices', 'Index for entry in %d-th dimension exceeds shape', i)
                end
                
                this.entries = containers.Map(subs2str(listIndices), listEntries(:,end));
                this.d       = d;
                this.shape   = shape;
            
            % N-dim matrix
            elseif nargin == 1 && isnumeric(varargin{1})
                denseEntries = varargin{1};
                this.shape = size(denseEntries);
                this.d     = length(this.shape);
                for i = 1:prod(this.shape)
                    subs = ind2subs(this.shape, i);
                    val = denseEntries(i);
                    this.addEntry(subs, val);
                end
                
            else
                error('Ndsparse:Constructor:InputArgsNotRecognized', 'Input arguments format not recognized. See documentation for usage.')
                
            end
            
            this.removeZeros();
            
        end
        
        function disp(this)
            % Pretty-print Ndsparse when called w/o semicolon or otherwise displayed
            rep = [num2str(this.d) '-d sparse tensor with ' num2str(this.nnz()) ' nonzero entries\n'];
            inds = keys(this.entries);
            for ind = inds
                rep = [rep '(' ind{1} ')\t' subs2str(this.entries(ind{1})) '\n'];
            end
            fprintf(rep)
        end
        
        function val = getEntry(this, pos)
            % Get entry at pos as vector of indices
            ind = subs2str(pos);
            if isKey(this.entries, ind)
                val = this.entries(ind);
            else
                val = 0;
            end
        end
        
        function addEntry(this, pos, val)
            % Add entry at pos as vector of indices, adding to element if already present
            ind = subs2str(pos);
            if isKey(this.entries, ind)
                this.entries(ind) = this.entries(ind) + val;
            else
                this.entries(ind) = val;
            end
        end
        
        function removeEntry(this, pos)
            % Remove entry at pos as vector of indices. Equivalent to setting a
            % pos equal to 0.
            ind = subs2str(pos);
            if isKey(this.entries, ind)
                this.entries.remove(ind);
            end
        end
        
        function numEntries = nnz(this)
            numEntries = this.entries.length;
        end
        
        function removeZeros(this)
            % Remove zeros entries, making it truly sparse
            indList = keys(this.entries);
            for i = indList
                ind = i{1};
                if this.entries(ind) == 0
                    this.entries.remove(ind);
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
            assert(isequal(sort(outDims), 1:length(outDims)), 'Ndsparse:ttt:InvalidOuterProdDimOrder', 'Outer product dimensions must contain integers 1:nLeftoverDims')
            
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
            indList1 = keys(this.entries);
            indList2 = keys(that.entries);
            
            % Perform contraction, only taking nonzero entries
            for i = indList1
                inds1 = i{1};
                pos1 = str2subs(inds1);
                inInds1 = pos1(in1);
                outInds1 = pos1(out1);
                
                for j = indList2
                    inds2 = j{1};
                    pos2 = str2subs(inds2);
                    inInds2 = pos2(in2);
                    outInds2 = pos2(out2);
                    
                    if all(inInds1 == inInds2)
                        outInds = [outInds1 outInds2];
                        pos = outInds(outDims);
                        val = this.entries(inds1) * that.entries(inds2);
                        prod.addEntry(pos, val);
                    end
                    
                end
            end
            
        end
        
        function transpose(this, permutation)
            % Transpose/permute dimensions according to permutation in-place.
            
            assert(length(permutation) == this.d, 'Ndsparse:transpose:InvalidPermutationDim', 'permutation vector doesn''t have the right number of dims')
            assert(isequal(sort(permutation), 1:this.d), 'Ndsparse:transpose:InvalidPermutationVector', 'permutation vector must contain the integer entries 1:d')
            
            listIndices = str2subs(keys(this.entries)');
            listIndices = listIndices(:,permutation); % perform permutation
            listIndexStrs = subs2str(listIndices);
            
            listValues = values(this.entries);
            
            this.entries = containers.Map(listIndexStrs, listValues);
            this.shape = this.shape(permutation);
        end
        
        function out = full(this)
            % Convert Ndsparse tensor to Matlab full N-dim matrix
            out = zeros(this.shape);
            indList = keys(this.entries);
            for i = indList
                ind = i{1};
                idx = subs2ind(this.shape, str2subs(ind));
                out(idx) = this.entries(ind); % linear indexing
            end
        end
        
        
    end
    
end