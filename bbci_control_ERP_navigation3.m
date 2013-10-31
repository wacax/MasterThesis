function [packet, state]= bbci_control_ERP_navigation(cfy_out, state, event, varargin)
%BBCI_CONTROL_ERP_SPELLER - Generate control signal for ERP-based Hex-o-Spell
%
%Synopsis:
%  PACKET= bbci_control_ERP_Speller(CFY_OUT, EVENT, SETTINGS, STATE)
%
%Arguments:
%  CFY_OUT - Output of the classifier
%  STATE - Internal state variable
%  EVENT - Structure that specifies the event (fields 'time' and 'desc')
%      that triggered the evaluation of this control.
%  SETTINGS - Structure specifying the following feedback paramters:
%    .nClasses - number of classes of stimuli from which to select
%    .nSequences - number of sequences (repetitions) that should be used
%        in order to determine one choice
%
%Output:
% PACKET: Variable/value list in a CELL defining the control signal that
%     is to be sent via UDP to the Speller application.
% STATE: Updated internal state variable
 
% 02-2013 Martijn Schreuder
 
opt = opt_proplistToStruct(varargin{:});
 
if isempty(state),
  state.forceMax = opt.forceMax;
  state.buffer = zeros(opt.nClasses, opt.bufferSize);
  state.counter = 0;
  state.medianValue = opt.medianValue;
  state.zeroValue = opt.zeroValue;
  if ~isempty(opt.window),
      window = opt.window(1,1:opt.bufferSize)/sum(opt.window(1,1:opt.bufferSize)); %normalize
      state.window = repmat(window, opt.nClasses,1);
  else
      state.window = ones(opt.nClasses, opt.bufferSize)/opt.bufferSize;
  end
  if ~isfield(opt, 'mapping') | isempty(opt.mapping),
      state.mapping = eye(opt.nClasses);
  else
      state.mapping = opt.mapping;
  end
  if ~isfield(opt, 'directions') | isempty(opt.directions),
    state.directions = ones(opt.nClasses,3);
  else
    state.directions = opt.directions;
  end
   
end
 
% this_cue= opt.mrk2feedback_fcn(event.desc);

cue = event.desc;
state.buffer(cue,:) = [state.buffer(cue, 2:end) cfy_out];
state.counter = state.counter + 1;

if ~mod(state.counter, opt.update_ival),
    %ver 1.-  
    %scores = mean(state.window .* state.buffer, 2)';
    %scores = median(state.window .* state.buffer, 2)';
    %scores
    %scores = sign(scores)./(1+exp(-state.zeroValue .* (abs(scores) - state.medianValue)));
    %scores
    
    %ver 5.- 
    scores = (state.window .* state.buffer);
    scores = sign(scores)./(1+exp(-state.zeroValue .* (abs(scores) - state.medianValue)));
    
    
    %ver 4.-
    %preScores = state.window .* state.buffer;
    %regularization = mean(preScores, 2);
    %[~, location] = max(kabs(preScores), [], 2);
    %for i = 1:size(preScores, 1)
    %  scores(i) = preScores(i , location(i));
    %end
    %scores = scores + regularization';
    %scores
    
    
    %ver 2.- 
    %scores = state.window .* state.buffer;
    %medianValue = 0.8;
    %zeroValue = 6;
    %scores = sum(sign(scores)./(1+exp(-zeroValue .* (abs(scores) - medianValue))), 2)';
    %scores = mean(sign(scores)./(1+exp(-zeroValue .* (abs(scores) - medianValue))), 2)';

    
    %ver 3.- 
    %scores = sum(state.window .* state.buffer, 2)';
    %scores = mean(state.window .* state.buffer, 2)';
    %range = abs(max(scores(:)) - min(scores(:)))/3
    %medianValue = 0.8;
    %zeroValue = 6;
    %scores = sign(scores)./(1+exp(-(state.zeroValue /range) .* (abs(scores) - state.medianValue*range)));
    
    %Packet prep with other mapping
    %ver 1.- Does not work
%     packet = zeros(1,size(state.mapping, 1));
%     for dim = 1:size(state.mapping, 1)
%         packet(dim) = mean([scores(state.mapping(dim, :) > 0) scores(state.mapping(dim, :) < 0)]);
%         if isempty(packet(dim)), packet(dim) = 0; end
%         %packet(dim) = mean(scores(state.mapping(dim, :) > 0));
%         %packet(dim) = min(scores(state.mapping(dim, :)));
%         if isnan(packet(dim))
%           packet(dim) = 0;
%         end
%     end
    
    %ver 2.- Not tested
    %packet = zeros(1,size(state.directions, 1));
    %for dim = 1:size(state.directions, 1),
    %    if isempty(scores),
    %        packet(dim) = 0;
    %    end
    %    packet(dim) = mean(scores(state.directions(dim, :) > 0));        
    %end
    
    %ver 3.-Current
    %for dim = 1:size(state.mapping,1),
    %    cl1 = min(scores(state.mapping(dim,:) > 0));
    %    cl2 = min(scores(state.mapping(dim,:) < 0));
    %    if isempty(cl1), cl1 = 0; end
    %    if isempty(cl2), cl2 = 0; end
    %    if abs(cl1-cl2) > .5,
    %      packet(dim) = min(1, max(cl1-cl2,-1));
    %    end
    %    packet
    %end
    
        %Packet prep
    %New Mapping
    packet = zeros(1,size(state.directions, 2));
%     for vecs = find(scores <= 0),
%       packet = packet + (abs(scores(vecs))*state.directions(vecs,:));
%     end
    [dum vecs] = min(scores);
    for i = 1:length(vecs),
      packet = packet + (abs(scores(vecs))*state.directions(vecs,:));
    end
    packet = packet / length(vecs);
    if state.forceMax,
      packet = sign(packet);
    end
    packet
   
   

    
    %Other Mapping
    %packet = zeros(1, size(state.directions, 2));
    %[dum vecs] = min(scores);
    %packet = packet + (abs(scores(vecs))*state.directions(vecs,:));
    %packet = []
    %packet 
    
    

else
    packet = [];
end