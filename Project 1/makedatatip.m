function varargout = makedatatip(hObj,index)
matlab_rel = version('-release');
MODIFY_MAKEDATATIP = ~(verLessThan('matlab', '8.4'));

% Check # of inputs
narginchk(2, 2)
nargoutchk(0, 1)

if length(hObj)~=1
  error('MAKEDATATIP:InvalidSize',...
    'HOBJ must be scalar.');
end

% Ensure hObj is valid target
if ~ishandle(hObj)
  error('MAKEDATATIP:InvalidHandle',...
    'HOBJ is an invalid handle object.');
end

isImage = strcmp(get(hObj, 'Type'), 'image'); %Determine if target is image

% Read data from hObj
try
  X = get(hObj,'XData');
  Y = get(hObj,'YData');
catch ME
  % Object must have an XData and YData property to be valid
  error('MAKEDATATIP:InvalidObjectType',...
    'Objects of class ''%s'' are not a valid targets for datatips.',...
    class(handle(hObj)))
end
try
  Z = get(hObj,'ZData');
catch ME
  % Many objects do not have a ZData property.  Some will work, some will
  % not.
  isImage = true;
end
% Ensure subscripts or indices are valid values and sizes
if isempty(index)
  return
elseif ~isnumeric(index)
  error('MAKEDATATIP:InvalidDataType',...
    'Subscript indices must be of numeric data type.')
elseif any(index(:)<1) ||...
    any(fix(index(:))~=index(:)) ||...
    any(isinf(index(:)))
  error('MAKEDATATIP:InvalidIndex',...
    'Subscript indices must be positive integers.')
elseif ~isvector(index) && ~any(size(index)==2)
  error('MAKEDATATIP:InvalidIndexMatrixSize',...
    'Subscript indices must be a vector or N-by-2 matrix.')
elseif (~isImage && isvector(X)) || size(index,2)~=2
  hDatatip = zeros(size(index));
  index = index(:);
  isLinear = true;
else
  hDatatip = zeros(size(index,1),1);
  isLinear = false;
end

% Get handle to datacursor mode object
hDataCursorMgr = datacursormode(ancestor(hObj,'figure'));

% Loop through each specified data point
for n = 1:size(index,1)
  
  % Create position vector
  if isImage && isLinear
    [i j] = ind2sub([X(2) Y(2)], index(n));
    pos = [i j 1];
  elseif isImage
    pos = [index(n, 1) index(n, 2) 1];
  elseif isempty(Z)
    if MODIFY_MAKEDATATIP
        pos = [X(index(n)) Y(index(n)) 0]; % For 2014b release or later
    else
        pos = [X(index(n)) Y(index(n))]; % Before 2014b release
    end
  elseif isLinear
    pos = [X(index(n)) Y(index(n)) Z(index(n))];
  else
    pos = [...
      X(index(n,1),index(n,2))...
      Y(index(n,1),index(n,2))...
      Z(index(n,1),index(n,2))];
  end
  
  % Create datatip
  hDatatip(n) = createDatatip(hDataCursorMgr, hObj);
  % Specify data cursor properties
    if ~MODIFY_MAKEDATATIP % Does not work in 2014b release or later
        if isImage
        set(get(hDatatip(n),'DataCursor'),'DataIndex',pos,...
          'TargetPoint',pos(1:2))
        else
        set(get(hDatatip(n),'DataCursor'),'DataIndex',index(n, :),...
          'TargetPoint',pos)
        end
    end
  
  % Specify datatip properties
  set(hDatatip(n),'Position',pos)
  
end

% Update all data cursors
updateDataCursors(hDataCursorMgr)

% Return handles if requested
if nargout==1
  varargout = {hDatatip};
end