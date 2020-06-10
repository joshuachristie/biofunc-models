% pTable should be a 2D table of the form:
% [s1 s2 s3 ... t p(s1,...t); % one s1,...,t combination
%  s1 s2 s3 ... t p(s1,...t); % another s1,...,t combination
%  s1 s2 s3 ... t p(s1,...t)]
% 
% All p's should add to 1

function normalised_PID_both_sources_present = redAsMinIComponentTotal(pTable)
    fprintf('\n');    
	numSources = size(pTable,2)-2;
	if (numSources <= 0)
		error('Not enough sources\n');
	end
	
	if (numSources > 2)
		error('Cannot handle > 2 sources');
	end
	
	pColumn = size(pTable,2);
	tColumn = size(pTable,2) - 1;
	if (abs(sum(pTable(:, pColumn)) - 1) > 0.00000000001)
		errorStr = sprintf('ps do not sum to 1 (%f)', sum(pTable(:,pColumn)));
		error(errorStr);
	end

	% Print headers:
	for s = 1 : numSources
		fprintf('s%d\t', s);
	end
	fprintf('t\tp\t');
	for s = 1 : numSources
		fprintf('p(s%d)\t', s);
	end
	fprintf('p(t)\t');
	for s = 1 : numSources
		fprintf('p(s%d|t)\t', s);
	end
	fprintf('p(S)\tp(S,t)\t');
	for s = 1 : numSources
		fprintf('i+(s%d)\ti-(s%d)\t', s, s);
	end
	fprintf('r+\t');
	for s = 1 : numSources
		fprintf('u+(s%d)\t', s);
	end
	fprintf('c+\t');
	fprintf('r-\t');
	for s = 1 : numSources
		fprintf('u-(s%d)\t', s);
	end
	fprintf('c-\t');
	fprintf('\n');
	rowsInPTable = size(pTable,1);
        
        % initialise vector to store information from s1=1, s2=1
        % (t={0,1})
        raw_PID_both_sources_present = zeros(2,8);
        
	for r = 1 : rowsInPTable
		% For each row, which is a unique sample configuration
		posInfo = zeros(1, numSources);
		negInfo = zeros(1, numSources);
		for s = 1 : numSources
			sourceVal = pTable(r,s);
			fprintf('%d\t', sourceVal);
		end
		targetVal = pTable(r,tColumn);
		fprintf('%d\t', targetVal);
		pTarget = sum(pTable(:,pColumn) .* (pTable(:,tColumn)==targetVal));
		fprintf('%.3f\t', pTable(r,pColumn));
		% Print pSource's
		pSources = zeros(1,numSources);
		iPosS = zeros(1,numSources);
		for s = 1 : numSources
			sourceVal = pTable(r,s);
			pSource = sum(pTable(:,pColumn) .* (pTable(:,s)==sourceVal));
			pSources(s) = pSource;
			iPosS(s) = - log2 (pSource);
			fprintf('%.3f\t', pSource);
		end
		fprintf('%.3f\t', pTarget);
		% Print pSourceGivenTarget's
		pSourceGivenTargets = zeros(1, numSources);
		iNegS = zeros(1,numSources);
		for s = 1 : numSources
			sourceVal = pTable(r,s);
			pSourceAndTarget = sum(pTable(:,pColumn) .* (pTable(:,s)==sourceVal).*(pTable(:,tColumn)==targetVal));
			% fprintf('%.3f\t', pSourceAndTarget);
			pSourceGivenTarget = pSourceAndTarget ./ pTarget;
			pSourcesGivenTarget(s) = pSourceGivenTarget;
			iNegS(s) = -log2 (pSourceGivenTarget);
			fprintf('%.3f\t', pSourceGivenTarget);
		end
		% Total information:
		pSourcesJoint = sum((sum(pTable(:,1:numSources) == repmat(pTable(r,1:numSources), rowsInPTable, 1), 2) == numSources).*pTable(:,pColumn));
		pSourcesJointAndTarget = sum((sum(pTable(:,1:(numSources+1)) == repmat(pTable(r,1:(numSources+1)), rowsInPTable, 1), 2) == (numSources+1)).*pTable(:,pColumn));
		pSourcesJointAndTargetVec(r) = pSourcesJointAndTarget;
        pSourcesJointGivenTarget = pSourcesJointAndTarget ./ pTarget;
		totalIPos = -log2(pSourcesJoint);
		totalINeg = -log2(pSourcesJointGivenTarget);
		fprintf('%.3f\t%.3f\t', pSourcesJoint, pSourcesJointAndTarget);
		% Print i+ and i- for each source
		for s = 1 : numSources
			fprintf('%.3f\t%.3f\t', iPosS(s), iNegS(s));
		end
		% Now do positive components
		rPos = min(iPosS);
        rPosVec(r) = rPos;
		fprintf('%.3f\t', rPos);
		uPoss = iPosS - rPos;
        uPosVec(:,r) = uPoss;
		for s = 1 : numSources
			fprintf('%.3f\t', uPoss(s))
		end
		% Just write the total higher order synegy terms:
		% higherOrderPos = totalIPos - max(iPosS); This isn't right
		% For 2 only!
		higherOrderPos = totalIPos - iPosS(1) - uPoss(2);
        higherOrderPosVec(r) = higherOrderPos;
		fprintf('%.3f\t', higherOrderPos);
		% Now do negative components
		rNeg = min(iNegS);
        rNegVec(r) = rNeg;
        fprintf('%.3f\t', rNeg);
		uNegs = iNegS - rNeg;
        uNegVec(:,r) = uNegs;
		for s = 1 : numSources
			fprintf('%.3f\t', uNegs(s))
        end
		% Just write the total higher order synegy terms:
		% higherOrderNeg = totalINeg - max(iNegS); This isn't right
		% For 2 only!
		higherOrderNeg = totalINeg - iNegS(1) - uNegs(2);
        higherOrderNegVec(r) = higherOrderNeg;
		fprintf('%.3f\n', higherOrderNeg);
        higherOrderComb = higherOrderPos - higherOrderNeg;
        
        % record terms if s1 = 1 and s2 = 1 (t = {0,1})
        if r < 3
            % [r+, u+(s1), u+(s2), c+, r-, u-(s1), u-(s2), c-]
            raw_PID_both_sources_present(r, :) = [rPos, uPoss(1), ...
                                uPoss(2), higherOrderPos, rNeg, ...
                                uNegs(1), uNegs(2), higherOrderNeg];
        end
        
    end

    rPosVec(isnan(rPosVec)) = 0;
    rPosVec(isinf(rPosVec)) = 0;
    rNegVec(isnan(rNegVec)) = 0;
    rNegVec(isinf(rNegVec)) = 0;

    uPosVec(isnan(uPosVec)) = 0;
    uPosVec(isinf(uPosVec)) = 0;
    uNegVec(isnan(uNegVec)) = 0;
    uNegVec(isinf(uNegVec)) = 0;

    higherOrderPosVec(isnan(higherOrderPosVec)) = 0;
    higherOrderPosVec(isinf(higherOrderPosVec)) = 0;
    higherOrderNegVec(isnan(higherOrderNegVec)) = 0;
    higherOrderNegVec(isinf(higherOrderNegVec)) = 0;
    
    RMinPos = sum(pSourcesJointAndTargetVec.*rPosVec);
    RMinNeg = sum(pSourcesJointAndTargetVec.*rNegVec);

    U1Pos = sum(pSourcesJointAndTargetVec.*uPosVec(1,:));
    U1Neg = sum(pSourcesJointAndTargetVec.*uNegVec(1,:));
    
    U2Pos = sum(pSourcesJointAndTargetVec.*uPosVec(2,:));
    U2Neg = sum(pSourcesJointAndTargetVec.*uNegVec(2,:));    
    
    CPos = sum(pSourcesJointAndTargetVec.*higherOrderPosVec);
    CNeg = sum(pSourcesJointAndTargetVec.*higherOrderNegVec);
    
    R = RMinPos - RMinNeg;
    U1 = U1Pos - U1Neg;
    U2 = U2Pos - U2Neg;
    C = CPos - CNeg;
    
    fprintf('\n');
    fprintf('U1\t');
    fprintf('U2\t');
    fprintf('R\t');
    fprintf('C\t\n');
    fprintf('%.3f\t', U1);
    fprintf('%.3f\t', U2);
    fprintf('%.3f\t', R);
    fprintf('%.3f\t\n\n', C);
    
    % calculate weighted decomposition for s1 = 1, s2 = 1, t =
    % {0,1} (i.e. top two lines of the decomposition)
    % order is R, U1, U2, C
    weighted_PID_both_sources_present = zeros(1,4);
    normalised_probs = [pTable(1,4), pTable(2,4)] / (pTable(1,4) + pTable(2,4));
    for i = 1:4
        weighted_PID_both_sources_present(i) = (raw_PID_both_sources_present(1, i) - ...
            raw_PID_both_sources_present(1, i + 4)) * normalised_probs(1) + ...
            (raw_PID_both_sources_present(2, i) - raw_PID_both_sources_present(2, i + 4)) ...
            * normalised_probs(2);
    end
    % normalise so that we can talk about decomposing the function metric
    % value (as opposed to talking about the amount of information in bits)
    normalised_PID_both_sources_present = ...
        weighted_PID_both_sources_present / sum(weighted_PID_both_sources_present);
end

