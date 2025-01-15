document.addEventListener('DOMContentLoaded', () => {
    function initializeTabs(containerSelector, tabSelector, contentSelector, defaultActive = true) {
        const container = document.querySelector(containerSelector);
        if (!container) {
            console.log(`Container not found: ${containerSelector}`);
            return;
        }

        const tabs = container.querySelectorAll(tabSelector);
        const tabContents = container.querySelectorAll(contentSelector);

        function activateTab(tabName) {
            tabContents.forEach(tab => tab.classList.remove('active'));
            tabs.forEach(tab => tab.classList.remove('active'));

            const activeContent = container.querySelector(`#${tabName}`);
            const activeTab = container.querySelector(`${tabSelector}[data-tab="${tabName}"]`);

            if (activeContent) {
                activeContent.classList.add('active');
            } else {
                console.log(`Content not found for tab: ${tabName}`);
            }

            if (activeTab) {
                activeTab.classList.add('active');
            } else {
                console.log(`Tab not found: ${tabName}`);
            }
        }

        tabs.forEach(tab => {
            tab.addEventListener('click', () => activateTab(tab.dataset.tab));
        });

        if (defaultActive && tabs.length > 0) {
            activateTab(tabs[0].dataset.tab);
        }
    }
    // Initialize BUSCO tabs
    initializeTabs('.busco-file', '.tab-link', '.tab');

    // Initialize tabs for QC files
    initializeTabs('.qc-file', '.tab-link', '.tab');

    // Initialize tabs for QUAST files
    initializeTabs('.quast-file', '.tab-link', '.tab');

    // Initialize annotation summary tabs
    initializeTabs('.annotation-summary', '.tab-link[data-tab^="summary-100-plus-"]', '.tab[id^="summary-100-plus-"]');
    initializeTabs('.annotation-summary', '.tab-link[data-tab^="summary-plus-"]', '.tab[id^="summary-plus-"]');

    // Show more/less functionality
    document.querySelectorAll('.show-more-btn').forEach(button => {
        button.addEventListener('click', () => {
            const section = document.getElementById(button.dataset.section);
            if (section) {
                section.classList.toggle('active');
                button.textContent = section.classList.contains('active') ? 'Show Less' : 'Show All';
            } else {
                console.log(`Section not found: ${button.dataset.section}`);
            }
        });
    });

    // Initialize main tabs for RSS MEME files
    const mainTabs = document.querySelectorAll('.rss-meme-file > .tabs .tab-link');
    const mainTabContents = document.querySelectorAll('.rss-meme-file > .tab');

    function activateMainTab(tabName) {
        mainTabContents.forEach(tab => tab.classList.remove('active'));
        mainTabs.forEach(tab => tab.classList.remove('active'));
        document.getElementById(tabName).classList.add('active');
        document.querySelector(`.rss-meme-file .tab-link[data-tab="${tabName}"]`).classList.add('active');
    }

    mainTabs.forEach(tab => {
        tab.addEventListener('click', function () {
            activateMainTab(this.dataset.tab);
        });
    });

    if (mainTabs.length > 0) {
        activateMainTab(mainTabs[0].dataset.tab);
    }

    // Initialize sub-tabs within each main tab of RSS MEME files
    mainTabContents.forEach(mainTab => {
        const subTabs = mainTab.querySelectorAll('.tabs .tab-link');
        const subTabContents = mainTab.querySelectorAll('.tab');

        function activateSubTab(tabName) {
            subTabContents.forEach(tab => {
                tab.classList.remove('active');
            });
            subTabs.forEach(tab => {
                tab.classList.remove('active');
            });
            document.getElementById(tabName).classList.add('active');
            document.querySelector(`.rss-meme-file .tab-link[data-tab="${tabName}"]`).classList.add('active');
        }

        subTabs.forEach(tab => {
            tab.addEventListener('click', function () {
                activateSubTab(this.dataset.tab);
            });
        });

        if (subTabs.length > 0) {
            activateSubTab(subTabs[0].dataset.tab);
        }
    });

    // Initialize table sorter
    $(document).ready(function () {
        $("table").tablesorter({
            theme: 'blue',
            widgets: ['zebra', 'columns'],
            widgetOptions: {
                zebra: ["even", "odd"],
                columns: ["primary", "secondary", "tertiary"]
            }
        });
    });
});