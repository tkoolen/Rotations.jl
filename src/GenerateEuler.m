function GenerateEuler()
%function GenerateEuler()
% Generate classes for various Euler angle combinations automatically

base = 'R';  % what to call the rotation matrix

% tait ordering scheme
tait_orders = cellstr(perms('XYZ'));  % extrinsic ordering
ntait = numel(tait_orders);

% proper euler ordering scheme
proper_orders = cellfun(@(order)([order(1:2), order(1)]), tait_orders, 'uniformoutput', false); % extrinsic ordering
nproper = numel(proper_orders);

% set up structures for tait order an proper euler ordering
tait_params = struct('schema',      'tait-byran',  ...
                     'type_name',   'EulerAngles', ... 
                     'order_type',  'TaitByranOrder', ...
                     'fields',      {{'theta_x', 'theta_y', 'theta_z'}}, ...
                     'ext_order',   {tait_orders}, ...
                     'Rext',        {cell(1,ntait)}, ...
                     'theta_sol',   {cell(1,ntait)});
                 
proper_params = struct('schema',      'proper-euler',  ...
                       'type_name',   'ProperEulerAngles', ... 
                       'order_type',  'ProperEulerOrder', ...
                       'fields',      {{'theta_1', 'theta_2', 'theta_3'}}, ...
                       'ext_order',   {proper_orders}, ...
                       'Rext',        {cell(1,nproper)}, ...
                       'theta_sol',   {cell(1,nproper)});    
                  
     
% generate symbolic rotation matrix
Rs = SymR(base);

% solve tait byron angles
tait_params   = SolveOrdering(tait_params, Rs);
proper_params = SolveOrdering(proper_params, Rs);

% now build an output file with the type info
fprintf('Creating: %s\n', fullfile(pwd, 'euler_types.jl'));
fid = fopen('euler_types.jl', 'w');
fprintf(fid, '%s', GenericHeader(''));
fprintf(fid, '\n\nusing FixedSizeArrays');
fprintf(fid, '\n\n%s', BuildOrderTypes(tait_params, proper_params));
fprintf(fid, '\n\n%s',  BuildTaitEulerType(tait_params));
fprintf(fid, '\n\n%s',  BuildProperEulerType(proper_params));
fid = fclose(fid);

% now build an output file with the transformation functions
fprintf('Creating: %s\n', fullfile(pwd, 'euler_conversions.jl'));
fid = fopen('euler_conversions.jl', 'w');
note = sprintf('# Working for Matrix -> EulerAngles follows the (hopefully) numerical stable method described in: %s\n', 'https://d3cw3dd2w32x2b.cloudfront.net/wp-content/uploads/2012/07/euler-angles1.pdf');
note = sprintf('%s# Note though that unlike the reference this package uses active (right handed) rotations\n', note);
fprintf(fid, '%s', GenericHeader(note));
fprintf(fid, '\n\n%s', BuildConversions(tait_params, base));
fprintf(fid, '\n\n%s', BuildConversions(proper_params, base));
fid = fclose(fid);
fprintf('Finished!\n');




function str = GenericHeader(note)
% function GenericHeader(note)

v = ver('symbolic');
str = sprintf('# File created: %s\n# using Matlab %s with %s version %s\n%s\n', datestr(now()), version('-release'), v.Name, v.Version, note);

function str = BuildConversions(params, base)
% function str = BuildConversions(params, base)
% function to build converions to / from Euler angles to rotation matices

hash = repmat('#', 1, 60);
str = sprintf('%s\n# Rotation matrix converions for Euler Angles (%s)\n%s\n', hash, params.schema, hash);

for i = 1:length(params.ext_order)
    
    % add a comment
    str = sprintf('%s\n\n#\n# Euler%s\n#\n', str, params.ext_order{i});
    
    % build untype to rotation matrix
    str = sprintf('%s\n%s', str, build_eulertorot(params, i, false));
    
    % build typed to rotation matrix
    str = sprintf('%s\n%s', str, build_eulertorot(params, i, true));
    
    % build untype from rotation matrix
    str = sprintf('%s\n\n%s', str, build_rottoeuler(params, i, false, base));
    
    % build typed from rotation matrix
    str = sprintf('%s\n%s', str, build_rottoeuler(params, i, true, base));
    
end

function str = build_rottoeuler(params, idx, eltyped, base)
% function str = build_eulertorot(params, idx, eltyped)
% create the function to turn a rotation matrix into a string

% typed or untyped output?
if (eltyped)
    str = sprintf('function rottoeuler{T}(::Type{%s{Euler%s,T}}, %s::RotMatrix)\n', params.type_name, params.ext_order{idx}, base);
    type_wrap = @(T, x)(sprintf('%s(%s)', T, x));  
    out_type = sprintf('%s{Euler%s,T}', params.type_name, params.ext_order{idx});
else
    str = sprintf('rottoeuler{T}(::Type{%s{Euler%s}}, %s::RotMatrix{T})', params.type_name, params.ext_order{idx}, base);
    str = sprintf('%s = rottoeuler(%s{Euler%s, T}, %s)\n', str, params.type_name, params.ext_order{idx}, base);
    return;
end

% replace R%i_%i style with R[i,j], also expand squares
theta_sol = cell(1, numel(fieldnames(params.theta_sol{idx})));
fnames = fieldnames(params.theta_sol{idx}); 
for i = 1:numel(theta_sol)
    theta_sol{i} = char(params.theta_sol{idx}.(fnames{i}));
    theta_sol{i} = regexprep(theta_sol{i}, [base,'(?<row>\d+)_(?<col>\d+)'], [base,'[$1, $2]']);
    theta_sol{i} = regexprep(theta_sol{i}, [base,'\[(?<row>\d+)\s*,\s*(?<col>\d+)\]\s*\^\s*2'], [base,'[$1, $2] * ', base, '[$1, $2]']);
end

% add the theta one working
str = sprintf('%s\n\tt1 = %s\n', str, theta_sol{1});
str = sprintf('%s\tct1, st1 = cos(t1), sin(t1)\n\n', str);
theta_sol{1} = 't1';

% type wrap solutions
for i = 1:numel(theta_sol)
    theta_sol{i} = type_wrap('T', theta_sol{i});
end

% put the solutiions in storage order
[~, idx] = ismember(fnames, params.fields);
theta_sol(idx) = theta_sol;

% adn write it
str = sprintf('%s\t%s(\n\t\t%s,\n\t\t%s,\n\t\t%s\n\t\t)\n', str, out_type, theta_sol{:});
str = sprintf('%send\n', str);

function str = build_eulertorot(params, idx, eltyped)
% function str = build_eulertorot(params, idx, eltyped)
% create the function to turn a rotation matrix into a string

base = 'ea';
% typed or untyped output?
if (eltyped)
    str = sprintf('function eulertorot{T,U}(::Type{RotMatrix{T}}, %s::%s{Euler%s,U})\n\n', base, params.type_name, params.ext_order{idx});
    type_wrap = @(T, x)(sprintf('%s(%s)', T, x));  
    str = sprintf('%s\twT = promote_type(T, U)', str);
else
    
    % allow spcifiny the Rotmatrix output without the element type
    str = sprintf('eulertorot{T}(::Type{RotMatrix}, %s::%s{Euler%s,T}) = eulertorot(%s::%s{Euler%s,T})\n\n', base, params.type_name, params.ext_order{idx}, base, params.type_name, params.ext_order{idx});
    
    % create a version that needs the first input
    str = sprintf('%sfunction eulertorot{T}(%s::%s{Euler%s,T})\n\n', str, base, params.type_name, params.ext_order{idx});
    type_wrap = @(T, x)(sprintf('%s', x));
end

% shortcuts
R = params.Rext{idx};
str = sprintf('%s\n', str);
for i = 1:numel(params.fields)
    fvar = type_wrap('wT', sprintf('%s.%s', base, params.fields{i}));
    str = sprintf('%s\tct%i, st%i = cos(%s), sin(%s)\n', str, i, i, fvar, fvar);
    
    % make substitions
    for j = 1:numel(params.Rext{idx})
        R(j) = subs(R(j), {cos(sym(params.fields{i})), sin(sym(params.fields{i}))}, {sprintf('ct%i', i), sprintf('st%i', i)});
    end
end

% and build it
str = sprintf('%s\n%s', str, mattofsa(R, @(x)type_wrap('T', x)));
str = sprintf('%s\nend\n', str);



function str = mattofsa(matrix, type_wrap)
% function str = mattofsa(R, type_wrap)
% build an fsa construction of R_str

% pad R to all the same size
matrix_str = cell(size(matrix));
for i = 1:numel(matrix)
    if iscell(matrix)
        el = char(matrix{i});
    else
        el = char(matrix(i));
    end
    
    % term cant start with minus then a space
    el = regexprep(el, '^\s*-\s*', '-');
    matrix_str{i} = type_wrap(el);
  
end
nchars = cellfun(@numel, matrix_str);
c_width = max(nchars, [], 1);

str = '';
for i = 1:size(matrix,1)
    if (i == 1)
        str = sprintf('%s\t@fsa([ ', str);
    else
        str = sprintf('%s\t       ', str);
    end
    str = sprintf('%s%s%s', str, matrix_str{i,1}, repmat(' ', 1, c_width(1) - numel(matrix_str{i,1})));
    for j = 2:size(matrix,2)
        str = sprintf('%s  %s', str, matrix_str{i,j}, repmat(' ', 1, c_width(j) - numel(matrix_str{i,j})));
    end
    if (i < size(matrix,1))
        str = sprintf('%s;\n', str);
    else
        str = sprintf('%s])\n', str);
    end
end
    



function header = BuildOrderTypes(tait_params, proper_params)
% function BuildOrderTypes(tait_order_name, tait_order, tait_int, euler_order_name, proper_order, proper_int)
% function to create the header for euler orders

% comment header
hash = repmat('#', 1, 30);
header = sprintf('%s\n# Euler Representations\n%s\n', hash, hash);

% abstract type for all headers
header = sprintf('%s\n\n"""\nEulerOrder specifies the extrinsic (axes of rotation are fixed) order of rotations\n"""\n', header);
header = sprintf('%sabstract EulerOrder     # Defines the order of rotations\n', header);

% loop over them
cparams = [tait_params, proper_params];
for i = 1:numel(cparams)
    header = sprintf('%s\n\n#\n# %s ordering\n#\n', header, cparams(i).schema);
    header = sprintf('%sabstract %s   <: EulerOrder\n', header, cparams(i).order_type);
    for j = 1:length(cparams(i).ext_order)
        header = sprintf('%s\n\n', header);
        header = sprintf('%s%s', header, OrderType(cparams(i), j));
    end
end

function header = OrderType(params, idx)
% function header = OrderType(ext_order, int_order, super_type, istait)
% Print the type for a specific order

if strncmpi(params.schema, 'tait', 4)
    field_ord = params.fields(axis_IDs(params.ext_order{idx}));
else
    field_ord = params.fields;
end
    
% the doc
header = sprintf('"""\n');
header = sprintf('%sEuler%s corresponds to:\n', header, params.ext_order{idx});
header = sprintf('%s\n\tExtrinsic: R = R%s(%s) * R%s(%s) * R%s(%s)\n', header, lower(params.ext_order{idx}(1)), field_ord{1}, lower(params.ext_order{idx}(2)), field_ord{2}, lower(params.ext_order{idx}(3)), field_ord{3});

% add the intrinsic form
header = sprintf('%s\nor equivalently,\n', header);
header = sprintf('%s\n\tIntrinsic: R = R%s''''(%s) * R%s''(%s) * R%s(%s)', header, lower(params.ext_order{idx}(3)), field_ord{3}, lower(params.ext_order{idx}(2)), field_ord{2}, lower(params.ext_order{idx}(1)), field_ord{1});
header = sprintf('%s (where the axis of rotation inherits from previous rotations)\n', header);

header = sprintf('%s"""\n', header);
header = sprintf('%stype Euler%s <: %s end\n', header, params.ext_order{idx}, params.order_type);


function header = BuildTaitEulerType(params)
% function header = BuildTaitEulerType(type_name, order_name, fields)
% function to create the header for euler angles with tait byran ordering

% comment header
hash = repmat('#', 1, 30);
header = sprintf('%s\n# Euler (Tait Byran) Angles\n%s\n', hash, hash);

% abstract type for all headers
header = sprintf('%s\n\n"""\n%s{Order <: %s, T <: AbstractFloat} where order specifies the extrinsic (axes of rotation are fixed) order of rotations.\n', header, params.type_name, params.order_type);
header = sprintf('%sSee subtypes(Rotations.%s) for a list of possible orders\n', header, params.order_type);
header = sprintf('%s\nFields:\n', header);
header = sprintf('%s\t%s   - Angle to rotate about the x axis (radians)\n', header, params.fields{1});
header = sprintf('%s\t%s   - Angle to rotate about the y axis (radians)\n', header, params.fields{2});
header = sprintf('%s\t%s   - Angle to rotate about the z axis (radians)\n', header, params.fields{3});
header = sprintf('%s"""\n', header);
header = sprintf('%simmutable %s{Order <: %s, T <: AbstractFloat} <: FixedVectorNoTuple{3, T}\n', header, params.type_name, params.order_type);
header = sprintf('%s\t%s::T\n', header, params.fields{1});
header = sprintf('%s\t%s::T\n', header, params.fields{2});
header = sprintf('%s\t%s::T\n', header, params.fields{3});
header = sprintf('%send\n', header);


function header = BuildProperEulerType(params)
% function header = BuildProperEulerType(params)
% function to create the header for euler angles with proper Euler ordering

% comment header
hash = repmat('#', 1, 30);
header = sprintf('%s\n# Proper Euler Angles\n%s\n', hash, hash);

% abstract type for all headers
header = sprintf('%s\n\n"""\n%s{Order <: %s, T <: AbstractFloat} where order specifies the extrinsic (axes of rotation are fixed) order of rotations.\n', header, params.type_name, params.order_type);
header = sprintf('%sSee subtypes(Rotations.%s) for a list of possible orders\n', header, params.order_type);
header = sprintf('%s\nFields:\n', header);
header = sprintf('%s\t%s   - angle to rotate about the first axis in the order (radians)\n', header, params.fields{1});
header = sprintf('%s\t%s   - angle to rotate about the second axis in the order (radians)\n', header, params.fields{2});
header = sprintf('%s\t%s   - angle to rotate about the final axis in the order (radians)\n', header, params.fields{3});
header = sprintf('%s"""\n', header);
header = sprintf('%simmutable %s{Order <: %s, T <: AbstractFloat} <: FixedVectorNoTuple{3, T}\n', header, params.type_name, params.order_type);
header = sprintf('%s\t%s::T\n', header, params.fields{1});
header = sprintf('%s\t%s::T\n', header, params.fields{2});
header = sprintf('%s\t%s::T\n', header, params.fields{3});
header = sprintf('%send\n', header);










function params = SolveOrdering(params, Rs)
% function params = SolveOrdering(params, Rs)
% solve the to and from rotation matrix problem

% do tait angles

for i = 1:numel(params.ext_order)
    disp(params.ext_order{i});
    [params.Rext{i}, R12, R3, sym_order] = GenEulerMatrix(params.schema, params.ext_order{i}, params.fields);
    params.theta_sol{i} = InverseRot(params.Rext{i}, R12, R3, Rs, sym_order);
    for j = 1:length(params.fields)
        fprintf('%s = %s\n', params.fields{j}, char(params.theta_sol{i}.(params.fields{j})));
    end
    disp(' ');
end



function res = InverseRot(Rf, R12, R3, Rs, params)
% function InverseRot(Rf, R12, R3, Rs, params)
%
% Rf     - symbolic rotation matrix for all rotations, Rf = R12 * R3
% R12    - symbolic rotation matrix after the first two rotations
% R3     - symbolic rotation matrix for the last rotation
% Rs     - symbolic result
% params - the symbols in Rf

% function to go back to angle from the rotation matrix

% build solution structure
var_str = {char(params(1)), char(params(2)), char(params(3))};
res = struct(var_str{1}, [], var_str{2}, [], var_str{3}, []); 

% find the frequency of cos and sin terms in Rf
freq = freq_count(Rf, params);

% This solution method should be more numerical stable than a symbolic solution

% first term - it should have a c * sin and a c * cos elements
ct_idx = find((sum(freq,2) == 2) & (freq(:,1) == 1));
st_idx = find((sum(freq,2) == 2) & (freq(:,2) == 1));
res.(var_str{1}) = atan_sol(Rf(st_idx) / Rf(ct_idx), Rs(st_idx) / Rs(ct_idx));

% second term term - start with the pure term
p_idx = find(sum(freq,2) == 1);
iscos = (freq(p_idx, 3) == 1);
ot_idx = find((sum(freq(:,3:end),2) == 2) & (freq(:,3+iscos) == 1));  % could use its interaction with theta1 as well
ot = sqrt(Rs(ot_idx(1))*Rs(ot_idx(1)) + Rs(ot_idx(2))*Rs(ot_idx(2)));
if iscos
    ds = Rf(p_idx) / cos(params(2));
    res.(var_str{2}) = atan2(ot, ds * Rs(p_idx));
else
    ns = Rf(p_idx) / sin(params(2));
    res.(var_str{2}) = atan2(ns * Rs(p_idx), ot);
end


% now solve for the third (possibly gimble locked) rotations
% for numerical stability, solve the system:
% Rf = R12 *  R3 = Rs ->
% R3 = R12' * Rs

% find the frequency of cos and sin terms in R12' and R3
R12d = R12';  
[freq12, R12sub] = freq_count(R12d, params(1:2));

% look for the row in R12 with only two terms
nt = reshape(sum(freq12,2), 3, 3);
ridx = find(sum(nt,2) == 2);

% find the non zero entries on that row in R3
[freq3, R3sub] = freq_count(R3, params(3), 2);
freq3 = freq3(ridx:3:end, :);  % remove rows that arent ridx
ct_idx = find(freq3(:,1) == 1);
st_idx = find(freq3(:,2) == 1);

% and solve 
sexpr = simple(solve(expand(R3sub(ridx, st_idx) - R12sub(ridx,:) * Rs(:, st_idx)), symvar(R3sub(ridx, st_idx)))); 
cexpr = simple(solve(expand(R3sub(ridx, ct_idx) - R12sub(ridx,:) * Rs(:, ct_idx)), symvar(R3sub(ridx, ct_idx)))); 
res.(var_str{3}) = atan2(sexpr, cexpr);



function [freq, R] = freq_count(R, params, skipped)
% function [freq, R] = freq_count(R, params, skipped)
% examing the contents of R
% this goes cos(params(1)), sin(params(1)), cos(params(2)) etc

if (nargin < 3)
    skipped = 0;
end

freq = zeros(9, 2*numel(params));
for i = 1:numel(params)
    R = subs(R, cos(params(i)), sprintf('ct%i',i+skipped));
    R = subs(R, sin(params(i)), sprintf('st%i',i+skipped));
end

% and count occurances in each term
for i = 1:numel(R)
    term_vars = symvar(R(i));
    for j = 1:numel(term_vars)
        str = char(term_vars(j));
        issin = str(1) == 's';
        var_id = str2double(str(end));
        freq(i, 2*(var_id-1-skipped)+1+issin) = 1;
    end
end




function sol = atan_sol(lhs, rhs)
% function to create the solution for varaints of  sin(x) / cos(x) = a / b
[nl, dl] = numden(lhs);
[nr, dr] = numden(rhs);
var = symvar(nl);

% need to invert both sides so its sin / cos?
if (isequal(dl, sin(var)) || isequal(dl, -sin(var)))
    [nl, dl, nr, dr] = deal(dl, nl, dr, nr);
end
nl_s = nl / sin(var);  % should be 1,-1
dl_s = dl / cos(var);  % should be 1,-1
if isequal(tan(var), simple(expand((nl_s*nl) / (dl_s*dl))))  % sanity check
    sol = simple(expand(atan2(nl_s*nr, dl_s*dr)));
else
    error('unknown solution format: %s', char(lhs))
end



function [Rext, R12, R3, params] = GenEulerMatrix(type, order, fields)
% function [Rext, R12, R3, params] = GenEulerMatrix(type, order, fields)
% generate the rotation mtrix from variable names

% frames
Rext = eye(3);
I = eye(3);

% generate symbolics for the fields
vec_expr = 'params = [';
for i = 1:length(fields)
    eval(sprintf('syms %s real', fields{i}))
    vec_expr = sprintf('%s%s, ', vec_expr, fields{i});
end
vec_expr(end-1:end) = '];';
eval(vec_expr);
    
% for tait byran update the field order to how they're applied
if strncmpi(type, 'tait', 4)
    params = params(axis_IDs(order));
end

% add assumptions (not helping)
assume((params >= -pi) & (params <= pi))

% axis indices
axes = axis_IDs(order);

% and create
for j = 3:-1:1
    R = ArbRot(params(j), I(:, axes(j)));
    Rext = R * Rext;
end

% also return the working parts
R12 = ArbRot(params(1), I(:, axes(1))) *  ArbRot(params(2), I(:, axes(2)));
R3 = ArbRot(params(3), I(:, axes(3)));


function R = ArbRot(theta, u)
% function R = ArbRot(theta, u)
% a rotation about u by theta
uc = [0,     -u(3),    u(2);
      u(3),  0,       -u(1);
     -u(2),  u(1),     0];
 
uo = u(:) * u(:)'; 
R = cos(theta) * eye(numel(u)) + sin(theta) * uc + (1 - cos(theta)) * uo;


function idx = axis_IDs(order)
% axis_IDs(order)
% returns 1-3 for 'XYZ'    
idx = uint8(order) - uint8('X') + 1;

function Rs = SymR(base)
% function Rs = SymR(base)
% build a symbolic rotation and try to add all known properties
Rs = sym(base, [3,3]);
for i = 1:numel(Rs)
    assume(Rs(i), 'real')
    assumeAlso((Rs(i) >= -1) & (Rs(i) <= 1))
end
assumeAlso(Rs * Rs' == eye(3))
assumeAlso(det(Rs) == 1)




% @doc """
% EulerZXY corresponds to:
% 
% 	Extrinsic: R = R(z) * R(x) * R(y)      
% 
% 	or equivilently,
% 
%     Intrinsic: R = R(y'') * R(x') * R(z)   (where axes of rotation inherit previous rotations)
% """ ->
% type EulerZXY <: EulerOrder end
