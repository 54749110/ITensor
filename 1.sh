echo "=== 网络连接诊断 ==="
echo "1. 测试22端口(SSH):"
timeout 3 nc -zv github.com 22 2>&1 | grep -E "(succeeded|timed out|refused|No route)"

echo -e "\n2. 测试443端口(HTTPS备用):"
timeout 3 nc -zv github.com 443 2>&1 | grep -E "(succeeded|timed out|refused|No route)"

echo -e "\n3. 测试SSH详细连接(前5行):"
timeout 5 ssh -vvvT git@github.com 2>&1 | grep -E "(Connecting to|Connection established|kex_exchange_identification|read: Connection|Permission denied)" | head -10

echo -e "\n4. 检查代理设置:"
echo "   http_proxy: ${http_proxy:-未设置}"
echo "   https_proxy: ${https_proxy:-未设置}"
echo "   all_proxy: ${all_proxy:-未设置}"
echo "   git代理: $(git config --global --get http.proxy 2>/dev/null || echo 未设置)"