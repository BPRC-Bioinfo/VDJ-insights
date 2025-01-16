document.addEventListener('DOMContentLoaded', function () {
    const sidebarContainer = document.getElementById('sidebar-container');

    if (sidebarContainer) {
        // Preload the sidebar CSS
        const link = document.createElement('link');
        link.rel = 'stylesheet';
        link.href = '../css/sidebar.css';
        link.onload = () => {
            fetchSidebarContent();
        };
        document.head.appendChild(link);
    }

    function fetchSidebarContent() {
        fetch('sidebar.html')
            .then(response => response.text())
            .then(data => {
                sidebarContainer.innerHTML = data;
                initializeSidebar();
            })
            .catch(error => console.error('Error loading sidebar:', error));
    }

    function initializeSidebar() {
        const toggleBtn = document.getElementById('toggle-btn');
        const sidebar = document.getElementById('sidebar');
        const mainContent = document.getElementById('main-content');
        const sidebarTitle = sidebar.querySelector('.sidebar-header h2'); // Select the h2 element inside sidebar-header
        const menuTexts = sidebar.querySelectorAll('.menu-text');
        const logo = sidebar.querySelector('.logo');

        function toggleSidebar(event) {
            event.preventDefault(); // Prevent any default action like form submission
            sidebar.classList.toggle('minimized');
            mainContent.classList.toggle('shifted');

            const isMinimized = sidebar.classList.contains('minimized');
            sidebarTitle.style.visibility = isMinimized ? 'hidden' : 'visible'; // Toggle visibility
            menuTexts.forEach(text => text.style.visibility = isMinimized ? 'hidden' : 'visible');
            logo.style.display = isMinimized ? 'none' : 'block'; // Hide/show logo
        }

        // Add event listener for toggle button
        if (toggleBtn) {
            toggleBtn.addEventListener('click', debounce(toggleSidebar, 100));
        }
    }

    // Debounce function to limit the rate of function execution
    function debounce(func, wait) {
        let timeout;
        return function (...args) {
            const later = () => {
                clearTimeout(timeout);
                func.apply(this, args);
            };
            clearTimeout(timeout);
            timeout = setTimeout(later, wait);
        };
    }
});
