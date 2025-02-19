% function index = computeIndex(trigger)
% [pos, typ] = ismember(trigger, [10 20]); 
% %pos is logical array, 1 if matches trigger values provided; 
% % typ indicates position of trigger element in provided trigger list, e.g. 1 2 3 for 9 13 10 respectively, or 0 if no match
% index.pos = find(pos); %returns indicies of elements that are true (=1)
% index.typ = typ(pos == 1) > 1; %select elements where pos == 1 and typ > 1, meaning any trigger 20; output is logical array 
% end


% 1 = no distractor, 2 = distractor
% 0 = midline target, 1 = lateral target
% last digit is distractor position
% set_size == 4:
%     mid_pos = [1,3]
%     lat_pos = [2,4] (2 right, 4 left)
% elif set_size == 6:
%     mid_pos = [1,4]
%     lat_pos = [2,3,5,6] (2,3 right; 5,6 left)
% elif set_size == 10:
%     mid_pos = [1,6]
%     lat_pos = [2,3,4,5,7,8,9,10] (2,3,4,5 right; 7,8,9,10 left)
% set_size=4
function index = computeIndex(trigger, set_size)
    if set_size == 4
        [pos, typ] = ismember(trigger, [102 104 100 110]); % 10: no-distractor, 20: distractor
        index.pos = find(pos); % Get positions of all distractor or no-distractor triggers
        typ_matched = typ(pos);
        index.typ = zeros(size(typ_matched));
        index.typ(typ_matched == 1) = 1; %dright
        index.typ(typ_matched == 2) = 2; %dleft
        index.typ(typ_matched > 2) = 0; %dnone
    elseif set_size == 6
        [pos, typ] = ismember(trigger, [202 203 205 206 100]); %[202 203 205 206 100]
        index.pos = find(pos); 
        typ_matched = typ(pos);
        index.typ = zeros(size(typ_matched));
        index.typ(typ_matched == 1 | typ_matched == 2 | typ_matched == 3 | typ_matched == 4) = 1;
        % index.typ(typ_matched == 1 | typ_matched == 2) = 1;
        %index.typ(typ_matched == 3 | typ_matched == 4) = 2;
        index.typ(typ_matched == 5) = 0;
    elseif set_size == 10
        [pos, typ] = ismember(trigger, [203 204 208 209 100]); %[202 203 204 205 207 208 209 211 100]
        index.pos = find(pos); 
        typ_matched = typ(pos);
        index.typ = zeros(size(typ_matched));
        index.typ(typ_matched >= 1 & typ_matched <= 2) = 1;
        index.typ(typ_matched >= 3 & typ_matched <= 4) = 2;
        index.typ(typ_matched == 5) = 0;
    end
end
% 
% function index = computeIndex(trigger)
% 
%     [pos, typ] = ismember(trigger, [202 204]); % 10: no-distractor, 20: distractor
%     index.pos = find(pos); % Get positions of all distractor or no-distractor triggers
%     index.typ = typ(pos == 1) ; % Logical array where type is left distractor
% end