function CommandLineInterface(args)
    if isempty(args)
        usage();
    end

    [opts args] = parse_args(struct(), args);

    if length(args) < 1
        err("Missing mode");
    end
    mode = args{1};
    args = args(2:end);
    
    switch mode
      case 'pipeline'
        opts = struct();
        opts.('min_domain_size') = '3008';
        opts.('min_n_domain_size') = '50000';
        opts.('win_size') = '32';
        
        [opts args] = parse_args(opts, args);
        if length(args) ~= 2
            err('Expecting input dir and output dir');
        end
        
        RunIsoPlotterPipeline(args{1}, ...
                              args{2}, ...
                              opt_num(opts, 'min_domain_size'), ...
                              opt_num(opts, 'min_n_domain_size'), ...
                              opt_num(opts, 'win_size'));
      case 'plotseglen'
        opts = struct();
        [opts args] = parse_args(opts, args);

        if length(args) ~= 2
            err('Expecting input file and output file');
        end
        
        PlotSegLengthGC(args{1}, args{2});

      case 'plotspatial'
        opts = struct();
        [opts args] = parse_args(opts, args);

        if length(args) ~= 3
            err('Expecting input file, chr #, and output file');
        end
        
        PlotSpatialGC(args{1}, ...
                      arg_num('chr #', args{2}), ...
                      args{3});

      otherwise
        err(['Invalid mode: ' mode])
    end
end

function [numval] = arg_num(argname, argval)
    numval = str2num(argval);
    if isempty(numval)
        err(['Invalid numeric value for argument ' argname ' (' argval ')']);
    end
end

function [numval] = opt_num(opts, optsym)
    numval = str2num(opts.(optsym));
    if isempty(numval)
        err(['Invalid numeric value for option ' optsym2name(optsym) ...
            ' (' opts.(optsym) ')']);
    end
end

function [sym] = optname2sym(optname)
    sym = strrep(optname, '-', '_');
end

function [name] = optsym2name(optsym)
    name = strrep(optsym, '_', '-');
end

function [opts args] = parse_args(options, args)
    n = length(args);
    i = 1;
    opts = options;
    while i <= n
        if strncmp('--', args{i}, 2)
            optname = args{i++}(3:end);
            optsym = optname2sym(optname);
            if i > n
                err(["Expecting argument for " optname]);
            end
            if ~isfield(opts, optsym)
                err(["Invalid option: " optname]);
            end
            optval = args{i++};
            opts.(optsym) = optval;
        elseif args{i}(1) == '-'
            err("Must use '--' for options");
        else
            break;
        end
    end
    
    args = args(i:end);
end

function err(msg)
    disp(msg)
    exit(1)
end

function usage()
    disp('usage: isoplotter pipeline [opts...] <input_dir> <output_dir>');
    disp('');
    disp('         --min-domain-size <int>');
    disp('         --min-n-domain-size <int>');
    disp('         --win-size <int>');
    disp('');
    disp('       isoplotter plotseglen <input_file> <output_file>');
    disp('');
    disp('       isoplotter plotspatial <input_file> <chr #> <output_file>');

    exit(0)
end

CommandLineInterface(argv())