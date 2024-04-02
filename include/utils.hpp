#include <map>

using namespace std;

template<typename K, typename V>
bool mapsAreEqual(const std::map<K, V>& map1, const std::map<K, V>& map2) {
	return map1.size() == map2.size() && std::equal(map1.begin(), map1.end(), map2.begin());
}

template<class K, class V> void printMap(map<K, V> m) {
    typename map<K, V>::iterator it = m.begin();
 
    while (it != m.end()) {
        cout << "Key: " << it->first << ", Value: " << it->second << endl;
        ++it;
    }
}
