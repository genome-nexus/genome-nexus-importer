#!/bin/bash
# Custom entrypoint wrapper to fix MongoDB Docker race condition
# The official entrypoint fails because the port isn't released fast enough between init shutdown and final startup.

set -e

# Check if this is a fresh database (init will run)
if [ ! -f "/data/db/.mongodb_initialized" ]; then
    echo "=== Custom Entrypoint: Fresh database detected ==="
    
    # Temporarily modify the original entrypoint behavior by setting a trap
    # Run init separately, then add delay, then start final mongod
    
    # First, run the original entrypoint which will:
    # 1. Start temp mongod
    # 2. Run init scripts
    # 3. Kill temp mongod
    # 4. Try to start final mongod (this is where race condition happens)
        
    # Run mongod for init, then restart with delay
    
    # Start mongod in background for initialization
    echo "=== Starting MongoDB for initialization ==="
    mongod --bind_ip 127.0.0.1 --fork --logpath /var/log/mongodb/init.log
    
    # Wait for mongod to be ready
    sleep 2
    
    # Run init scripts manually
    echo "=== Running initialization scripts ==="
    for f in /docker-entrypoint-initdb.d/*; do
        case "$f" in
            *.sh)
                echo "Running $f"
                . "$f"
                ;;
            *.js)
                echo "Running $f"
                mongosh --quiet localhost:27017 "$f"
                ;;
        esac
    done
    
    # Mark as initialized
    touch /data/db/.mongodb_initialized
    
    # Stop init mongod gracefully
    echo "=== Stopping initialization MongoDB ==="
    mongod --shutdown || true
    
    # Wait for port to be fully released
    echo "=== Waiting for port release (5 seconds) ==="
    sleep 5
    
    # Verify port is free
    while netstat -tuln 2>/dev/null | grep -q ":27017 " || ss -tuln 2>/dev/null | grep -q ":27017 "; do
        echo "Port 27017 still in use, waiting..."
        sleep 1
    done
    
    echo "=== Starting final MongoDB ==="
fi

# Start the final mongod normally
exec mongod --bind_ip_all
