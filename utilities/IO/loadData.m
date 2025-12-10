function data = loadData(filename)
    loaded = load(filename, 'data');
    data = loaded.data;
end