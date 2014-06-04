function newMap = copyMap(oldMap)

keyType = oldMap.KeyType;
valueType = oldMap.ValueType;

newMap = containers.Map('KeyType',keyType,'ValueType',valueType);

oldKeys = keys(oldMap);
for i = 1:length(oldKeys)
    key = oldKeys{i};
    val = oldMap(key);
    newMap(key) = val;
end