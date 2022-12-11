#include <bits/stdc++.h>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
 
using namespace std;
using namespace __gnu_pbds;

template <typename K, typename V, typename Comp = std::less<K>>
using ordered_map = tree<K, V, Comp, rb_tree_tag, tree_order_statistics_node_update>;
template <typename K, typename Comp = std::less<K>>
using ordered_set = ordered_map<K, null_type, Comp>;

template <typename K, typename V, typename Comp = std::less_equal<K>>
using ordered_multimap = tree<K, V, Comp, rb_tree_tag, tree_order_statistics_node_update>;
template <typename K, typename Comp = std::less_equal<K>>
using ordered_multiset = ordered_multimap<K, null_type, Comp>;

#define ll long long
#define nl "\n"
#define mod 1'000'000'007
#define cin(vec) for (auto &i : vec) cin >> i
#define cout(vec) for (auto &i : vec) cout << i << " ";cout << '\n';
#define all(v) v.begin(), v.end()
#define rall(v) v.rbegin(),v.rend()
#define fixed(n) fixed << setprecision(n)
#define ceil(n, m) (((n) / (m)) + ((n) % (m) ? 1 : 0))
#define Time cerr << "Time Taken: " << (float)clock()/CLOCKS_PER_SEC << " secs\n";

void HamokshaIsinTheHouse()
{
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    #ifndef ONLINE_JUDGE
    freopen("in.txt", "r", stdin);
    freopen("out.txt", "w", stdout);
    #endif
}

void solve()
{       
    
}

int main()
{
    HamokshaIsinTheHouse();
    int t = 1;
    // cin >> t;  
    while(t--)
        solve();
    Time 
    return 0;
}