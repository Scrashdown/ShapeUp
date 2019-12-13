//=============================================================================
//
//   General mesh viewer using NanoGUI and OpenGL.
//   Supports mesh loading, orbit viewing, highlighting vertices,
//   grabbing and dragging mesh vertices with the mouse and selecting
//   vertices.
//
//   "Digital 3D Geometry Processing"
//
//   Gaspard Zoss, Christopher Brandt
//
//   Copyright (C) 2019 by Computer Graphics and Geometry Laboratory,
//         EPFL
//
//-----------------------------------------------------------------------------
#pragma once

#include <nanogui/opengl.h>
#include <nanogui/glutil.h>
#include <nanogui/screen.h>
#include <nanogui/window.h>
#include <nanogui/layout.h>
#include <nanogui/popupbutton.h>
#include <nanogui/label.h>
#include <nanogui/button.h>
#include <nanogui/textbox.h>
#include <nanogui/tabwidget.h>
#include <surface_mesh/Surface_mesh.h>
#include <random> // Added just to test temperature there

#if defined(__GNUC__)
#  pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#endif
#if defined(_WIN32)
#  pragma warning(push)
#  pragma warning(disable: 4457 4456 4005 4312)
#endif

#if defined(_WIN32)
#  pragma warning(pop)
#endif
#if defined(_WIN32)
#  if defined(APIENTRY)
#    undef APIENTRY
#  endif
#  include <windows.h>
#  define M_PI 3.1415926535897932384626433832795028841971
#endif

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;
using std::pair;
using std::to_string;
using std::min;
using std::max;
using namespace surface_mesh;
using surface_mesh::Color;
using namespace nanogui;

typedef unsigned long Index;

constexpr float INITIAL_GRAB_RADIUS = 20;
constexpr float DEFAULT_SCALE = 0.5;

enum VertexStatus {
    Default, Selected, Pinned, Movable
};

class Viewer : public nanogui::Screen {
public:
    Viewer(std::string title, bool (*pre_draw_callback)(Viewer*) = nullptr, bool (*mesh_load_callback)(Viewer*) = nullptr, bool (*pre_mesh_load_callback)(Viewer*) = nullptr)
        :
        nanogui::Screen(Eigen::Vector2i(1024, 768), title) {

        initGUI();
        initShaders();

		m_pre_draw_callback = pre_draw_callback;
		m_mesh_load_callback = mesh_load_callback;
        m_pre_mesh_load_callback = pre_mesh_load_callback;
    }

    void loadMesh(string filename) {
        if (m_pre_mesh_load_callback) {
            if (!m_pre_mesh_load_callback(this)) {
                std::cout << "Error on callback before loading mesh!" << std::endl;
                return;
            }
        }

        if (!mesh.read(filename)) {
            std::cerr << "Mesh not found, exiting." << std::endl;
            exit(-1);
        }

        cout << "Mesh "<< filename << " loaded." << endl;
        n_vertices = mesh.n_vertices();
        cout << "# of vertices : " << n_vertices << endl;
        n_faces = mesh.n_faces();
        cout << "# of faces : " << n_faces << endl;
        n_edges = mesh.n_edges();
        cout << "# of edges : " << n_edges << endl;

        /** NEW: build vertex lookup table **/
        m_vertex_lookup_table.resize(mesh.n_vertices());
        for (const auto v : mesh.vertices()) {
            m_vertex_lookup_table[v.idx()] = v;
        }

        mesh_center = computeCenter(&mesh);
        float dist_max = 0.0f;
        for (auto v: mesh.vertices()) {
            if ( distance(mesh_center, mesh.position(v))> dist_max) {
                dist_max = distance(mesh_center, mesh.position(v));
            }
        }

        std::cout << "Scale before: " << dist_max;

        // Rescale mesh, such that dist_max = SCALE;
        for (auto v : mesh.vertices()) {
            mesh.position(v) = (mesh.position(v) - mesh_center) * (DEFAULT_SCALE / dist_max) + mesh_center;
        }
        dist_max = DEFAULT_SCALE;

        float dist_after = 0.0f;
        for (auto v : mesh.vertices()) {
            if (distance(mesh_center, mesh.position(v)) > dist_after) {
                dist_after = distance(mesh_center, mesh.position(v));
            }
        }

        std::cout << ", scale after: " << dist_after << std::endl; 

        mCamera.arcball = Arcball(2.);
        mCamera.arcball.setSize(mSize);
        mCamera.modelZoom = 2/dist_max;
        mCamera.modelTranslation = -Vector3f(mesh_center.x, mesh_center.y, mesh_center.z);

		meshProcess();

        if (m_mesh_load_callback) {
            if (!m_mesh_load_callback(this)) {
                std::cout << "Error on callback after loading mesh!" << std::endl;
            }
        }
	}

	void updateShaderVertices(const MatrixXf& vPos, bool forced = false) {
        if ((!m_reupload_vertices) || forced) {
            m_updated_shader_verts = vPos;
            m_reupload_vertices = true;
        }
	}

    void updateShaderNormals(const MatrixXf& vNormals, bool forced = false) {
        if ((!m_reupload_normals) || forced) {
            m_updated_shader_normals = vNormals;
            m_reupload_normals = true;
        }
    }

    void updateVertexStatus(const Eigen::Matrix<int, 1, -1>& vStatus) {
        m_current_vertex_status = vStatus;
        updateVertexStatusVisualization();
    }

    const Eigen::Matrix<int, 1, -1>& getVertexStatus() {
        return m_current_vertex_status;
    }

    /**NEW! Update color matrix for temperature visualization */
    void setTemperatureColor() {
        Surface_mesh::Vertex_property<Scalar> v_temperature =
                mesh.vertex_property<Scalar>("v:temperature", 0.0);
        Surface_mesh::Vertex_property<surface_mesh::Color> v_color_temperature =
                mesh.vertex_property<surface_mesh::Color>("v:color_temperature",
                                                          surface_mesh::Color(1.0f, 1.0f, 1.0f));
        for (auto v : mesh.vertices()) {

            set_color(v, value_to_color(v_temperature[v], c_min_value, c_max_value), v_color_temperature);
            m_color_temperature_.col(v.idx()) << v_color_temperature[v].x,
                    v_color_temperature[v].y,
                    v_color_temperature[v].z;
        }

        m_reupload_vertex_colors = true;

    }

        /** NEW: Updates temperature of selected vertices with provided temperature**/
    /**
     * Sets the temperature of the selected vertices to the float passed as argument
     * @param temperature Float value with which to update temperature of selected vertices
     */
    void setSelectedVerticesTemperature(const float temperature){
        Surface_mesh::Vertex_property<Scalar> v_temperature =
                mesh.vertex_property<Scalar>("v:temperature",0.0);
        Surface_mesh::Vertex_property<surface_mesh::Color> v_color_temperature =
                mesh.vertex_property<surface_mesh::Color>("v:color_temperature",
                                                          surface_mesh::Color(1.0f, 1.0f, 1.0f));

        /** The loop is only made on the selected vertices*/
        for(auto id : m_selectedVertices){
            const auto v = m_vertex_lookup_table[id];

            /** Update the relevant temperature quantities at last */
            v_temperature[v] = temperature;
            set_color(v, value_to_color(temperature, c_min_value, c_max_value), v_color_temperature);
            m_color_temperature_.col(id) <<  v_color_temperature[v].x,
                    v_color_temperature[v].y,
                    v_color_temperature[v].z;
            m_reupload_vertex_colors = true;
        }
    }

    /** NEW: Toggles the source status for selected vertex **/
    /**
     * Among the selected vertices, make all non-source vertices as vertices and make source as non-sources
     * Also updates the corresponding display nodes (sources are displayed as red dots, others are not displayed outside of selection)
     */
    void toggleTemperatureSourceSelectedVertices(){
        Surface_mesh::Vertex_property<bool> v_is_source = mesh.vertex_property<bool>("v:is_source", false);
        auto vtcs = mesh.vertices().begin();
        for(auto id: m_selectedVertices){
            const auto v = m_vertex_lookup_table[id];

            /** Flip the source boolean*/
            v_is_source[v] =! v_is_source[v];
            /** Updates the display: a source is displayed in red, otherwise it is not displayed*/
            if(v_is_source[v]){
                m_updated_vertex_sources.col(id) = Vector3f( 1.0 ,0,0);
                m_updated_vertex_selections.col(id) = Vector3f( 1.0 ,0,0);
            } else {
                m_updated_vertex_sources.col(id) = Vector3f(0.0, 0, 0);
                m_updated_vertex_selections.col(id) = Vector3f(0.0, 0, 0);
            }

            m_reupload_vertex_selections = true;
        }
    }

    void clearSources() {
        Surface_mesh::Vertex_property<bool> v_is_source = mesh.vertex_property<bool>("v:is_source", false);

        for(auto v: mesh.vertices()){
            v_is_source[v] = false;
        }
        m_reupload_vertex_selections = true;
    }

    void updateVertexStatusVisualization() {
        m_updated_vertex_selections.setZero(3, n_vertices);
        for (Index i = 0; i < std::max(m_updated_vertex_selections.cols(), m_updated_vertex_selections.cols()); i++) {
            switch (m_current_vertex_status(0, i)) {
            case VertexStatus::Movable:
                m_updated_vertex_selections.col(i) = Vector3f(0.133, 0.694, 0.298);
                break;
            case VertexStatus::Pinned:
                m_updated_vertex_selections.col(i) = Vector3f(0.5, 0.0, 0.0);
                break;
            case VertexStatus::Selected:
                m_updated_vertex_selections.col(i) = Vector3f(0.0, 0.635, 0.909);
                break;
            default:
                m_updated_vertex_selections.col(i) = m_updated_vertex_sources.col(i); /** The default display accounts for source display this way **/
            }
        }
        updateVertexSelectionVisualization();
        m_reupload_vertex_selections = true;
    }

    void updateVertexSelectionVisualization() {
        if (m_updated_vertex_selections.rows() != 3 || m_updated_vertex_selections.cols() != n_vertices) {
            m_updated_vertex_selections.setZero(3, n_vertices);
        }
        for (Index selVInd : m_selectedVertices) {
            m_updated_vertex_selections.col(selVInd) = Vector3f(0.0, 0.635, 0.909);
        }
        m_reupload_vertex_selections = true;
    }

    void setFloorHeight(float floorHeight) {
        m_floorHeight = floorHeight;
        m_floorHeightChanged = true;
    }

    void showFloor(bool show) {
        m_showFloor = show;
    }

    void meshProcess() {
        Point default_normal(0.0, 1.0, 0.0);
        Surface_mesh::Vertex_property<Point> vertex_normal =
                mesh.vertex_property<Point>("v:normal");
        mesh.update_face_normals();
        mesh.update_vertex_normals();

        /** NEW! Vertex temperature attribute **/
        Surface_mesh::Vertex_property<Scalar> v_temperature =
                mesh.vertex_property<Scalar>("v:temperature",0.0);

        Surface_mesh::Vertex_property<surface_mesh::Color> v_color_temperature =
                mesh.vertex_property<surface_mesh::Color>("v:color_temperature",
                                                          surface_mesh::Color(1.0f, 1.0f, 1.0f));

        /** Random fill of the v_temperature array**/
        std::default_random_engine generator;
        std::uniform_int_distribution<int> distribution(0,255);

        for(auto v: mesh.vertices()){
            v_temperature[v] = distribution(generator);
        }

        /** Converts temperature scalar to color vector **/
        color_coding(v_temperature, &mesh, v_color_temperature,255);

        int j = 0;
        MatrixXf mesh_points(3, n_vertices);
        MatrixXu indices(3, n_faces);
        m_color_temperature_ = MatrixXf(3, n_vertices);
        m_color_temperature_.setZero();

        for(auto f: mesh.faces()) {
            vector<float> vv(3.0f);
            int k = 0;
            for (auto v: mesh.vertices(f)) {
                vv[k] = v.idx();
                ++k;
            }
            indices.col(j) << vv[0], vv[1], vv[2];
            ++j;
        }

        // Create big matrices to send the data to the GPU with the required
        // format
        MatrixXf normals_attrib(3, n_vertices);

        j = 0;
        for (auto v: mesh.vertices()) {
            mesh_points.col(j) << mesh.position(v).x,
                                  mesh.position(v).y,
                                  mesh.position(v).z;

            normals_attrib.col(j) << vertex_normal[v].x,
                                     vertex_normal[v].y,
                                     vertex_normal[v].z;
            /** NEW: Temperature update**/
            m_color_temperature_.col(j) <<  v_color_temperature[v].x,
                                            v_color_temperature[v].y,
                                            v_color_temperature[v].z;

            ++j;
        }
        m_updated_shader_verts = mesh_points;

        MatrixXf vertex_color(3, n_vertices);
        vertex_color.setZero();
        for (int i = 0; i < n_vertices; i++) vertex_color.col(i) = Vector3f(0.98, 0.59, 0.04);

        mShader.bind();
        mShader.uploadIndices(indices);
        mShader.uploadAttrib("position", mesh_points);
        mShader.uploadAttrib("normal", normals_attrib);
        mShader.uploadAttrib("color", vertex_color);

        /** NEW: Temperature attrib upload **/
        mShader.uploadAttrib("temperature_color", m_color_temperature_);

        // Initialize floor geom and shader
        int floor_grid_length = 50;
        int num_floor_points = floor_grid_length * floor_grid_length;
        int num_floor_faces = (floor_grid_length - 1) * (floor_grid_length - 1) * 2;
        m_floorPoints.resize(3, num_floor_points);
        MatrixXf floor_normals(3, num_floor_points);
        MatrixXu floor_indices(3, num_floor_faces);
        float floor_min_xy = -100;
        float floor_max_xy = 100;
        float floor_cell_length = (float)(floor_max_xy - floor_min_xy) / (float)floor_grid_length;
        for (int x = 0; x < floor_grid_length; x++) {
            for (int y = 0; y < floor_grid_length; y++) {
                m_floorPoints(0, x * floor_grid_length + y) = floor_min_xy + floor_cell_length * x;
                m_floorPoints(1, x * floor_grid_length + y) = m_floorHeight;
                m_floorPoints(2, x * floor_grid_length + y) = floor_min_xy + floor_cell_length * y;
                floor_normals(0, x * floor_grid_length + y) = 0;
                floor_normals(1, x * floor_grid_length + y) = 1;
                floor_normals(2, x * floor_grid_length + y) = 0;
            }
        }
        for (int x = 0; x < floor_grid_length - 1; x++) {
            for (int y = 0; y < floor_grid_length - 1; y++) {
                int cell_ind = (x * (floor_grid_length - 1) + y) * 2;
                floor_indices(0, cell_ind) = x * floor_grid_length + y;
                floor_indices(1, cell_ind) = (x + 1) * floor_grid_length + y;
                floor_indices(2, cell_ind) = x * floor_grid_length + (y + 1);
                floor_indices(0, cell_ind + 1) = (x + 1) * floor_grid_length + y;
                floor_indices(1, cell_ind + 1) = (x + 1) * floor_grid_length + (y + 1);
                floor_indices(2, cell_ind + 1) = x * floor_grid_length + (y + 1);

            }
        }
        m_numFloorFaces = floor_indices.cols();
        mFloorShader.bind();
        mFloorShader.uploadIndices(floor_indices);
        mFloorShader.uploadAttrib("position", m_floorPoints);
        mFloorShader.uploadAttrib("normal", floor_normals);

        MatrixXf selected(3, n_vertices);
        selected.setZero();
        m_updated_vertex_selections.setZero(3, n_vertices);
        m_updated_vertex_sources.setZero(3, n_vertices);
        m_current_vertex_status.setZero(1, n_vertices);
        mSelectedVertexShader.bind();
        mSelectedVertexShader.uploadIndices(indices);
        mSelectedVertexShader.uploadAttrib("position", mesh_points);
        mSelectedVertexShader.uploadAttrib("selected", selected);

        MatrixXu simple_indices(3, 1);
        simple_indices << 0, 1, 2;
        MatrixXf zeros(3, 3);
        zeros.setZero();
        mSelectionQuadShader.bind();
        mSelectionQuadShader.uploadIndices(simple_indices);
        mSelectionQuadShader.uploadAttrib("position", zeros);
    }

    void initGUI() {
        window = new Window(this, "Controls");
        window->setPosition(Vector2i(15, 15));
        window->setLayout(new GroupLayout());

        PopupButton *popupBtn = new PopupButton(window, "Open a mesh", ENTYPO_ICON_EXPORT);
        Popup *popup = popupBtn->popup();
		popup->setLayout(new GroupLayout());

        Button* b = new Button(popup, "Eight");
        b->setCallback([this,popupBtn]() {
            loadMesh("../data/eight.off");
            popupBtn->setPushed(false);
        });

        b = new Button(popup, "Small Sphere");
        b->setCallback([this,popupBtn]() {
            loadMesh("../data/small_sphere.obj");
            popupBtn->setPushed(false);
        });

        b = new Button(popup, "Max-Planck");
        b->setCallback([this, popupBtn]() {
            loadMesh("../data/max.off");
            popupBtn->setPushed(false);
        });

        b = new Button(popup, "Open from disk");
        b->setCallback([this, popupBtn] {
            std::string load_file = file_dialog(
                { {"off", "OFF File"}, {"obj", "Obj File"} }, false);
            std::cout << load_file << std::endl;
            loadMesh(load_file);
            popupBtn->setPushed(false);
        });

        new Label(window, "Display Control", "sans-bold");

        b = new Button(window, "Wireframe");
        b->setFlags(Button::ToggleButton);
        b->setChangeCallback([this](bool wireframe) {
            this->wireframe =! this->wireframe;
        });

        /** NEW: Button for temperature **/
        b = new Button(window, "Temperature");
        b->setFlags(Button::ToggleButton);
        b->setChangeCallback([this](bool temperature){
            this->temperature = !this->temperature;
        });

        /** NEW: Toggle button for mesh surface display **/
        b = new Button(window, "Display Mesh Surface");
        b->setFlags(Button::ToggleButton);
        b->setPushed(displaySurface);
        b->setChangeCallback([this](bool displaySurface) {
            this->displaySurface = !this->displaySurface;
        });

        /** NEW: Box for temperature change **/
        Label* iterations_label = new Label(window, "Set Temperature: ");
        FloatBox<float>* iterations_box = new FloatBox<float>(window);
        iterations_box->setEditable(true);
        iterations_box->setCallback([this](float temperature) {
            setSelectedVerticesTemperature(temperature);
        });

        /** NEW: Button to do/undo source status of an edge **/
        b = new Button(window, "Make/unmake source");
        b->setFlags(Button::NormalButton);
        b->setCallback([this](void){
            toggleTemperatureSourceSelectedVertices();
        });

        /** NEW! Button for clearing all sources at once **/
        b = new Button(window, "Clear sources");
        b->setFlags(Button::NormalButton);
        b->setCallback([this](void){
            clearSelection();
            clearSources();
        });

        performLayout();
    }

    void initShaders() {
        // Shaders
        mShader.init(
            "a_simple_shader",

            /* Vertex shader */
            "#version 330\n"
            "uniform mat4 MV;\n"
            "uniform mat4 P;"
            "uniform int color_mode;\n"

            "in vec3 position;\n"
            "in vec3 normal;\n"
            "in vec3 color;\n"
            "in vec3 temperature_color;\n"

            "out vec3 fcolor;\n"
            "out vec3 fnormal;\n"
            "out vec3 view_dir;\n"
            "out vec3 light_dir;\n"

            "void main() {\n"
            "    vec4 vpoint_mv = MV * vec4(position, 1.0);\n"
            "    gl_Position = P * vpoint_mv;\n"
            "    if (color_mode == 1) {\n"
            "        fcolor = temperature_color;\n"
            "    } else {\n"
            "        fcolor = color;\n"
            "    }\n"
            "    fnormal = mat3(transpose(inverse(MV))) * normal;\n"
            "    light_dir = vec3(0.0, 3.0, 3.0) - vpoint_mv.xyz;\n"
            "    view_dir = -vpoint_mv.xyz;\n"
            "}",

            /* Fragment shader */
            "#version 330\n"
            "uniform vec3 intensity;\n"
            "uniform int color_mode;\n"

            "in vec3 fcolor;\n"
            "in vec3 fnormal;\n"
            "in vec3 view_dir;\n"
            "in vec3 light_dir;\n"

            "out vec4 color;\n"

            "void main() {\n"
            "    vec3 c = vec3(0.0);\n"
            "    if(color_mode == 0){ \n"
            "        c += vec3(1.0)*vec3(0.18, 0.1, 0.1);\n"
            "        vec3 n = normalize(fnormal);\n"
            "        vec3 v = normalize(view_dir);\n"
            "        vec3 l = normalize(light_dir);\n"
            "        float lambert = dot(n,l);\n"
            "        if(lambert > 0.0) {\n"
            "            c += vec3(1.0)*vec3(0.9, 0.5, 0.5)*lambert;\n"
            "            vec3 v = normalize(view_dir);\n"
            "            vec3 r = reflect(-l,n);\n"
            "            c += vec3(1.0)*vec3(0.8, 0.8, 0.8)*pow(max(dot(r,v), 0.0), 90.0);\n"
            "        }\n"
            "        c *= fcolor;\n"
            "    } else {\n"
            "        c = fcolor; \n"
            "    }\n"
            "    if (intensity == vec3(0.0)) {\n"
            "        c = intensity;\n"
            "    }\n"
            "    color = vec4(c, 1.0);\n"
            "}"
        );

        mFloorShader.init(
            "floor_shader",

            /* Vertex shader */
            "#version 330\n"
            "uniform mat4 MV;\n"
            "uniform mat4 P;\n"

            "in vec3 position;\n"
            "in vec3 normal;\n"

            "out vec3 fnormal;\n"
            "out vec3 view_dir;\n"
            "out vec3 light_dir;\n"

            "void main() {\n"
            "    vec4 vpoint_mv = MV * vec4(position, 1.0);\n"
            "    gl_Position = P * vpoint_mv;\n"
            "    fnormal = mat3(transpose(inverse(MV))) * normal;\n"
            "    light_dir = vec3(0.0, 3.0, 3.0) - vpoint_mv.xyz;\n"
            "    view_dir = -vpoint_mv.xyz;\n"
            "}",

            /* Fragment shader */
            "#version 330\n"
            "uniform vec3 intensity;\n"

            "in vec3 fnormal;\n"
            "in vec3 view_dir;\n"
            "in vec3 light_dir;\n"

            "out vec4 color;\n"

            "void main() {\n"
            "    vec3 c = vec3(0.0);\n"
            "    c += vec3(1.0)*vec3(0.18, 0.1, 0.1);\n"
            "    vec3 n = normalize(fnormal);\n"
            "    vec3 v = normalize(view_dir);\n"
            "    vec3 l = normalize(light_dir);\n"
            "    float lambert = dot(n,l);\n"
            "    if(lambert > 0.0) {\n"
            "        c += vec3(1.0)*vec3(0.9, 0.5, 0.5)*lambert;\n"
            "        vec3 v = normalize(view_dir);\n"
            "        vec3 r = reflect(-l,n);\n"
            "        c += vec3(1.0)*vec3(0.8, 0.8, 0.8)*pow(max(dot(r,v), 0.0), 90.0);\n"
            "    }\n"
            "    c *= vec3(0.23, 0.29, 0.4);\n"
            "    if (intensity == vec3(0.0)) {\n"
            "        c = intensity;\n"
            "    }\n"
            "    color = vec4(c, 1.0);\n"
            "}"
        );

        mSelectedVertexShader.init(
            "selected_vertex_shader",
            /* Vertex shader */
            "#version 330\n\n"

            "uniform mat4 MV;\n"
            "uniform mat4 P;\n"

            "in vec3 position;\n"
            "in vec3 selected;\n"

            "out vec3 vertexSelected;\n"

            "void main() {\n"
            "    vertexSelected = selected;\n"
            "    gl_Position = P * (MV * vec4(position, 1.0));\n"
            "}",

            /* Fragment shader */
            "#version 330\n\n"

            "out vec4 color;\n"

            "in vec3 diamondColor;\n"

            "void main() {\n"
            "    color = vec4(diamondColor, 1);\n"
            "}",


            /* Geometry shader */
            "#version 330\n\n"

            "layout (triangles) in;\n"
            "layout (line_strip, max_vertices = 15) out;\n"

            "uniform mat4 MV;\n"
            "uniform mat4 P;\n"

            "in vec3 vertexSelected[];\n"

            "out vec3 diamondColor;\n"

            "void creatediamond(int index) {\n"
            "   diamondColor = vertexSelected[index];\n"
            "   if (diamondColor.x + diamondColor.y + diamondColor.z < 0.00001) return;"
            "   gl_Position = gl_in[index].gl_Position;\n"
            "   gl_Position.x = gl_Position.x + 0.035f;\n"
            "   EmitVertex();\n"
            "   gl_Position.x = gl_Position.x - 0.035f;\n"
            "   gl_Position.y = gl_Position.y + 0.035f;\n"
            "   EmitVertex();\n"
            "   gl_Position.x = gl_Position.x - 0.035f;\n"
            "   gl_Position.y = gl_Position.y - 0.035f;\n"
            "   EmitVertex();\n"
            "   gl_Position.x = gl_Position.x + 0.035f;\n"
            "   gl_Position.y = gl_Position.y - 0.035f;\n"
            "   EmitVertex();\n"
            "   gl_Position.x = gl_Position.x + 0.035f;\n"
            "   gl_Position.y = gl_Position.y + 0.035f;\n"
            "   EmitVertex();\n"
            "   EndPrimitive();\n"
            "}\n"

            "void main() {\n"
            "   creatediamond(0);\n"
            "   creatediamond(1);\n"
            "   creatediamond(2);\n"
            "}"
        );

        mSelectionQuadShader.init(
            "selection_shader",
            /* Vertex shader */
            "#version 330\n\n"

            "in vec3 position;\n"

            "void main() {\n"
            "    gl_Position = vec4(position.x, position.y, 0. , 1.0);\n"
            "}",

            /* Fragment shader */
            "#version 330\n\n"

            "out vec4 color;\n"

            "void main() {\n"
            "    color = vec4(1, 1, 1, 1);\n"
            "}",


            /* Geometry shader */
            "#version 330\n\n"

            "layout (triangles) in;\n"
            "layout (line_strip, max_vertices = 5) out;\n"
            
            "void main() {\n"
            "   gl_Position = gl_in[0].gl_Position;\n"
            "   EmitVertex();\n"
            "   gl_Position = gl_in[1].gl_Position;\n"
            "   EmitVertex();\n"
            "   gl_Position = gl_in[2].gl_Position;\n"
            "   EmitVertex();\n"
            "   gl_Position.x = gl_in[2].gl_Position.x;\n"
            "   gl_Position.y = gl_in[0].gl_Position.y;\n"
            "   EmitVertex();\n"
            "   gl_Position = gl_in[0].gl_Position;\n"
            "   EmitVertex();\n"
            "   EndPrimitive();\n"
            "}"
        );
    }

    ~Viewer() {
        mShader.free();
        mSelectedVertexShader.free();
    }

    Point computeCenter(Surface_mesh *mesh) {
        Point center = Point(0.0f);

        for (auto v: mesh->vertices()) {
            center += mesh->position(v);
        }

        return center/mesh->n_vertices();
    }

    virtual bool keyboardEvent(int key, int scancode, int action, int modifiers) {
        if (Screen::keyboardEvent(key, scancode, action, modifiers)) {
            return true;
        }
        if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
            setVisible(false);
            return true;
        }
        return false;
    }

    virtual void draw(NVGcontext *ctx) {
        /* Draw the user interface */
        Screen::draw(ctx);
    }

    Vector2f getScreenCoord() {
        Vector2i pos = mousePos();
        return Vector2f(2.0f * (float)pos.x() / width() - 1.0f,
                        1.0f - 2.0f * (float)pos.y() / height());
    }

    void repaint() {
        //glfwPostEmptyEvent();
    }

    const std::vector<Index>& getSelectedVertices() {
        return m_selectedVertices;
    }

    void setSelectedVertices(const std::vector<Index>& vertInds) {
        m_selectedVertices = vertInds;
        updateVertexStatusVisualization();
    }


    virtual void drawContents() {
        using namespace nanogui;

		/* Pre-draw callback */
		if (m_pre_draw_callback) {
			if (!m_pre_draw_callback(this)) {
				return;
			}
		}

		/* Draw the window contents using OpenGL */
		mShader.bind();

		if (m_reupload_vertices) {
			mShader.uploadAttrib("position", m_updated_shader_verts);
		}

        if (m_reupload_normals) {
			mShader.uploadAttrib("normal", m_updated_shader_normals);
		}

        /** NEW: Temperature reupload**/
        if(m_reupload_vertex_colors){
            mShader.uploadAttrib("temperature_color", m_color_temperature_);
        }

        Eigen::Matrix4f model, view, proj;
        computeCameraMatrices(model, view, proj);

        Matrix4f mv = view*model;
        Matrix4f p = proj;

        /* MVP uniforms */
        mShader.setUniform("MV", mv);
        mShader.setUniform("P", p);

        /* Setup OpenGL (making sure the GUI doesn't disable these */
        glEnable(GL_DEPTH_TEST);
        glDisable(GL_CULL_FACE);

        /* Render everything */
        if (wireframe) {
            glEnable(GL_POLYGON_OFFSET_FILL);
            glPolygonOffset(1.0, 1.0);
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        }

        Vector3f colors(0.98, 0.59, 0.04);
        mShader.setUniform("intensity", colors);
        mShader.drawIndexed(GL_TRIANGLES, 0, n_faces);

        if (wireframe) {
            glDisable(GL_POLYGON_OFFSET_FILL);
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            colors << 0.0, 0.0, 0.0;
            mShader.setUniform("intensity", colors);
            mShader.drawIndexed(GL_TRIANGLES, 0, n_faces);
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        }

        /** NEW! Temperature update in the shader **/
        if(temperature){
            mShader.setUniform("color_mode", 1);
        } else {
            mShader.setUniform("color_mode", 0);
        }

        if (m_showFloor) {
            mFloorShader.bind();
            if (m_floorHeightChanged) {
                m_floorPoints.row(1).setConstant(m_floorHeight);
                mFloorShader.uploadAttrib("position", m_floorPoints);
                m_floorHeightChanged = false;
            }
            Vector3f colors(0.98, 0.59, 0.04);
            mFloorShader.setUniform("MV", mv);
            mFloorShader.setUniform("P", p);
            mFloorShader.setUniform("intensity", colors);
            mFloorShader.drawIndexed(GL_TRIANGLES, 0, m_numFloorFaces);
        }

        if (true) {
            mSelectedVertexShader.bind();
            if (m_reupload_vertex_selections) {
                mSelectedVertexShader.uploadAttrib("selected", m_updated_vertex_selections);
            }
            if (m_reupload_vertices) {
                mSelectedVertexShader.uploadAttrib("position", m_updated_shader_verts);
            }
            mSelectedVertexShader.setUniform("MV", mv);
            mSelectedVertexShader.setUniform("P", p);
            mSelectedVertexShader.drawIndexed(GL_TRIANGLES, 0, n_faces);
        }

        if (m_isSelecting) {
            mSelectionQuadShader.bind();
            float winWidth = mSize.x();
            float winHeight = mSize.y();
            float x1 = (m_selectionStart.x() / winWidth) * 2 - 1;
            float x2 = (m_selectionEnd.x() / winWidth) * 2 - 1;
            float y1 = -((m_selectionStart.y() / winHeight) * 2 - 1);
            float y2 = -((m_selectionEnd.y() / winHeight) * 2 - 1);
            MatrixXf curSelQuad(3, 3);
            curSelQuad <<
                x1, x1, x2,
                y1, y2, y2,
                0,  0,  0;
            mSelectionQuadShader.uploadAttrib("position", curSelQuad);
            mSelectionQuadShader.drawIndexed(GL_TRIANGLES, 0, 1);
        }

        m_reupload_normals = false;
        m_reupload_vertex_selections = false;
        m_reupload_vertices = false;
        /** NEW: reset of temperature boolean reupload flag**/
        m_reupload_vertex_colors = false;
    }

    bool scrollEvent(const Vector2i &p, const Vector2f &rel) {
        if (!Screen::scrollEvent(p, rel)) {
            mCamera.zoom = max(0.1, mCamera.zoom * (rel.y() > 0 ? 1.1 : 0.9));
            repaint();
        }
        return true;
    }

    void reset() {
        // reset all the components
        // recenter the mesh (maybe keep the original mesh somewhere so that if
        // we modify - smooth or else - it we can restore the original one.)
    }

    bool mouseMotionEvent(const Vector2i &p, const Vector2i &rel,
                          int button, int modifiers) {
        if (!Screen::mouseMotionEvent(p, rel, button, modifiers)) {
            if (m_isGrabbing) {
                grabMove(p);
            } else if (m_isSelecting) {
                selectionMove(p);
                repaint();
            } else if (mCamera.arcball.motion(p)) {
                repaint();
            } else if (mTranslate) {
                Eigen::Matrix4f model, view, proj;
                computeCameraMatrices(model, view, proj);
                float zval = nanogui::project(Vector3f(mesh_center.x,
                                                       mesh_center.y,
                                                       mesh_center.z),
                                              view * model, proj, mSize).z();
                Eigen::Vector3f pos1 = nanogui::unproject(
                        Eigen::Vector3f(p.x(), mSize.y() - p.y(), zval), view * model, proj, mSize);
                Eigen::Vector3f pos0 = nanogui::unproject(
                        Eigen::Vector3f(mTranslateStart.x(), mSize.y() -
                           mTranslateStart.y(), zval), view * model, proj, mSize);
                mCamera.modelTranslation = mCamera.modelTranslation_start + (pos1-pos0);
                repaint();
            }
        }
        return true;
    }

    bool mouseButtonEvent(const Vector2i &p, int button, bool down, int modifiers) {
        if (!Screen::mouseButtonEvent(p, button, down, modifiers)) {
            // In grab mode the left mouse click selects/drags vertices of the mesh,
            if (button == GLFW_MOUSE_BUTTON_1) {
                if (modifiers == GLFW_MOD_ALT) {
                    grabButton(p, down);
                } if (modifiers == GLFW_MOD_CONTROL) {
                    selectButton(p, down);
                } else if ( modifiers == 0) {
                    mCamera.arcball.button(p, down);
                }
            }
            if (button == GLFW_MOUSE_BUTTON_2 ||
                  (button == GLFW_MOUSE_BUTTON_1 && modifiers == GLFW_MOD_SHIFT)) {
                mCamera.modelTranslation_start = mCamera.modelTranslation;
                mTranslate = true;
                mTranslateStart = p;
            }
        }

        if (button == GLFW_MOUSE_BUTTON_1 && !down) {
            grabButton(p, down);
            selectButton(p, down);
            mCamera.arcball.button(p, down);
        }

        if (!down) {
            mDrag = false;
            mTranslate = false;
        }
        return true;
    }

    void setGrabCallbacks(void (*grabCallback)(const std::vector<Index>&,const std::vector<Vector3f>&), void (*grabReleaseCallback)()) {
        m_grab_callback = grabCallback;
        m_grab_release_callback = grabReleaseCallback;
    }

	Surface_mesh* getMesh() {
		return &mesh;
	}

    const std::vector<Surface_mesh::Vertex>& getVertexLookupTable() const { return m_vertex_lookup_table; }

    void clearSelection() {
        m_selectedVertices.clear();
        updateVertexStatusVisualization();
    }

    /** Method from the code framework for the lecture
     * "Digital 3D Geometry Processing"
     * Gaspard Zoss, Alexandru Ichim
     * Copyright (C) 2016 by Computer Graphics and Geometry Laboratory,
     * EPF Lausanne
     * @param prop The property array for which to associate colors
     * @param mesh The target mesh
     * @param color_prop The color property array to fill
     * @param bound The color bound
     */
    void color_coding(Surface_mesh::Vertex_property<Scalar> prop, Surface_mesh *mesh,
                      Surface_mesh::Vertex_property<surface_mesh::Color> color_prop, int bound=20) {
        // Get the value array
        std::vector<Scalar> values = prop.vector();

        // discard upper and lower bound
        unsigned int n = values.size()-1;
        unsigned int i = n / bound;
        std::sort(values.begin(), values.end());
        Scalar min_value = values[i], max_value = values[n-1-i];
        c_min_value = min_value;
        c_max_value = max_value;

        // map values to colors
        for (auto v: mesh->vertices())
        {
            set_color(v, value_to_color(prop[v], min_value, max_value), color_prop);
        }
    }

    /** Method from the code framework for the lecture
     * "Digital 3D Geometry Processing"
     * Gaspard Zoss, Alexandru Ichim
     * Copyright (C) 2016 by Computer Graphics and Geometry Laboratory,
     * EPF Lausanne
     * @param v Target vertex (id to update in the array)
     * @param col The color to update in the array (value to update)
     * @param color_prop Target array to update with color at target vertex's id
     */
    void set_color(Surface_mesh::Vertex v, const surface_mesh::Color& col,
                                   Surface_mesh::Vertex_property<surface_mesh::Color> color_prop)
    {
        color_prop[v] = col;
    }

    /** Method from the code framework for the lecture
     * "Digital 3D Geometry Processing"
     * Gaspard Zoss, Alexandru Ichim
     * Copyright (C) 2016 by Computer Graphics and Geometry Laboratory,
     * EPF Lausanne
     * @param value Scalar value to translate to an RGB triplet
     * @param min_value Lower bound of the color
     * @param max_value Upper bound of the color
     * @return The RGB-converted scalar value
     */
    surface_mesh::Color value_to_color(Scalar value, Scalar min_value, Scalar max_value) {
        Scalar v0, v1, v2, v3, v4;
        v0 = min_value + 0.0/4.0 * (max_value - min_value);
        v1 = min_value + 1.0/4.0 * (max_value - min_value);
        v2 = min_value + 2.0/4.0 * (max_value - min_value);
        v3 = min_value + 3.0/4.0 * (max_value - min_value);
        v4 = min_value + 4.0/4.0 * (max_value - min_value);

        surface_mesh::Color col(1.0f, 1.0f, 1.0f);

        if (value < v0) {
            col = surface_mesh::Color(0, 0, 1);
        } else if (value > v4) {
            col = surface_mesh::Color(1, 0, 0);
        } else if (value <= v2) {
            if (value <= v1) { // [v0, v1]
                Scalar u =  (value - v0) / (v1 - v0);
                col = surface_mesh::Color(0, u, 1);
            } else { // ]v1, v2]
                Scalar u = (value - v1) / (v2 - v1);
                col = surface_mesh::Color(0, 1, 1-u);
            }
        } else {
            if (value <= v3) { // ]v2, v3]
                Scalar u = (value - v2) / (v3 - v2);
                col = surface_mesh::Color(u, 1, 0);
            } else { // ]v3, v4]
                Scalar u = (value - v3) / (v4 - v3);
                col = surface_mesh::Color(1, 1-u, 0);
            }
        }
        return col;
    }



private:

    void grabButton(const Vector2i& mousePos, bool down) {
        if (m_isGrabbing && !down) {
            grabFree();
        }

        if (!m_isGrabbing && down) {
            grabSelect(mousePos);
        }
    }

    void grabFree() {
        m_isGrabbing = false;
        m_grabSelection.clear();
        if (m_grab_release_callback) m_grab_release_callback();
    }

    void grabSelect(const Vector2i& mousePos)
    {
    	if ( !m_currentMVP.hasNaN() && m_updated_shader_verts.rows() > 0 && m_grab_callback) {
    		int winWidth = mSize.x();
    		int winHeight = mSize.y();

    		m_grabSelection.clear();
    		m_grabVertPos.clear();
    		m_grabOrigPos.setZero();
    		m_grabOrigProj.setZero();
            // Check which projections of used vertices are close to mouse
    		for (Index i = 0; i < m_updated_shader_verts.cols(); i++) {
    			Vector3f vert = m_updated_shader_verts.col(i);
    			Vector4f vertProjection(vert(0), vert(1), vert(2), 1.f);
    			vertProjection = m_currentMVP * vertProjection;
    			Vector2i screenPos;
    			screenPos(0) = (((vertProjection(0) / vertProjection(3)) + 1.f) / 2.f) * (float)winWidth;
    			screenPos(1) = (((vertProjection(1) / -vertProjection(3)) + 1.f) / 2.f) * (float)winHeight;
    			float distToMouse = (screenPos - mousePos).norm();
    			if (distToMouse < m_grabRadius) {
    				m_grabSelection.push_back(i);
    				m_grabVertPos.push_back(vert);
    				m_grabOrigProj += vertProjection;
    			}
    		}
    		if (m_grabSelection.size() > 0) {
    			m_grabOrigProj *= 1.f / (float)m_grabSelection.size();
                m_grabOrigProj(0) = ((2.f * (float)mousePos.x() / (float)winWidth) - 1.f) * m_grabOrigProj(3);
                m_grabOrigProj(1) = ((2.f * (float)mousePos.y() / (float)winHeight) - 1.f) * (-m_grabOrigProj(3));
                m_grabOrigPos = (m_currentMVP.inverse() * m_grabOrigProj).template segment<3>(0);

                m_isGrabbing = true;
    			m_grab_callback(m_grabSelection, m_grabVertPos);
    		}
    		else {
    			m_isGrabbing = false;
    		}
    	}

    	m_lastMousePosition = mousePos;
    }

    void startSelecting(const Vector2i& mousePos) {
        m_isSelecting = true;
        m_selectionStart = mousePos;
        m_selectionEnd = mousePos;
    }

    void stopSelecting(const Vector2i& mousePos) {
        m_selectionEnd = mousePos;
        m_selectedVertices.clear();
        int winWidth = mSize.x();
        int winHeight = mSize.y();
        for (Index i = 0; i < m_updated_shader_verts.cols(); i++) {
            Vector3f vert = m_updated_shader_verts.col(i);
            Vector4f vertProjection(vert(0), vert(1), vert(2), 1.f);
            vertProjection = m_currentMVP * vertProjection;
            Vector2i screenPos;
            screenPos(0) = (((vertProjection(0) / vertProjection(3)) + 1.f) / 2.f) * (float)winWidth;
            screenPos(1) = (((vertProjection(1) / -vertProjection(3)) + 1.f) / 2.f) * (float)winHeight;
            if (screenPos(0) >= min(m_selectionStart.x(), m_selectionEnd.x())
                && screenPos(0) <= max(m_selectionStart.x(), m_selectionEnd.x())
                && screenPos(1) >= min(m_selectionStart.y(), m_selectionEnd.y())
                && screenPos(1) <= max(m_selectionStart.y(), m_selectionEnd.y())) {
                m_selectedVertices.push_back(i);
            }
        }
        updateVertexStatusVisualization();
        m_isSelecting = false;
    }

    void selectButton(const Vector2i& mousePos, bool down) {
        if (m_isSelecting && !down) {
            stopSelecting(mousePos);
        }

        if (!m_isSelecting && down) {
            startSelecting(mousePos);
        }
    }

    void selectionMove(const Vector2i& mousePos) {
        m_selectionEnd = mousePos;
    }

    void grabMove(const Vector2i& mousePos)
    {
        if ( m_isGrabbing && m_grabSelection.size() > 0 && !m_currentMVP.hasNaN() && m_updated_shader_verts.rows() > 0 && m_grab_callback) {
            int winWidth = mSize.x();
    		int winHeight = mSize.y();

    		Vector4f newProj = m_grabOrigProj;
    		newProj(0) = ((2.f * (float)mousePos.x() / (float)winWidth) - 1.f) * newProj(3);
    		newProj(1) = ((2.f * (float)mousePos.y() / (float)winHeight) - 1.f) * (-newProj(3));
    		Vector3f newPos = (m_currentMVP.inverse() * newProj).template segment<3>(0);

    		Vector3f trans = newPos - m_grabOrigPos;

    		std::vector<Vector3f> newVertPos;
    		for (Index v = 0; v < m_grabVertPos.size(); v++) {
    			newVertPos.push_back(m_grabVertPos[v] + trans);
    		}
    		m_grab_callback(m_grabSelection, newVertPos);
    	}

    	m_lastMousePosition = mousePos;
    }

    struct CameraParameters {
        nanogui::Arcball arcball;
        float zoom = 1.0f, viewAngle = 45.0f;
        float dnear = 0.05f, dfar = 100.0f;
        Eigen::Vector3f eye = Eigen::Vector3f(0.0f, 0.0f, 5.0f);
        Eigen::Vector3f center = Eigen::Vector3f(0.0f, 0.0f, 0.0f);
        Eigen::Vector3f up = Eigen::Vector3f(0.0f, 1.0f, 0.0f);
        Eigen::Vector3f modelTranslation = Eigen::Vector3f::Zero();
        Eigen::Vector3f modelTranslation_start = Eigen::Vector3f::Zero();
        float modelZoom = 1.0f;
    };

    CameraParameters mCamera;
    bool mTranslate = false;
    bool mDrag = false;
    Vector2i mTranslateStart = Vector2i(0,0);

    void computeCameraMatrices(Eigen::Matrix4f &model,
                               Eigen::Matrix4f &view,
                               Eigen::Matrix4f &proj) {

        view = nanogui::lookAt(mCamera.eye, mCamera.center, mCamera.up);

        float fH = std::tan(mCamera.viewAngle / 360.0f * M_PI) * mCamera.dnear;
        float fW = fH * (float) mSize.x() / (float) mSize.y();

        proj = nanogui::frustum(-fW, fW, -fH, fH, mCamera.dnear, mCamera.dfar);
        model = mCamera.arcball.matrix();
        model *= nanogui::scale(Eigen::Vector3f::Constant(mCamera.zoom * mCamera.modelZoom));
        model *= nanogui::translate(mCamera.modelTranslation);

        m_currentMVP = proj * view * model;
    }

    // Variables for the viewer
    nanogui::GLShader mShader;
    nanogui::GLShader mFloorShader;
    nanogui::GLShader mSelectedVertexShader;
    nanogui::GLShader mSelectionQuadShader;
    nanogui::Window *window;
    nanogui::Arcball arcball;

    Point mesh_center = Point(0.0f, 0.0f, 0.0f);
    Surface_mesh mesh;

    /** NEW: tightly packed array for finding a vertex given its index **/
    std::vector<Surface_mesh::Vertex> m_vertex_lookup_table;

    enum COLOR_MODE : int { NORMAL = 0, VALENCE = 1, CURVATURE = 2 };
    enum CURVATURE_TYPE : int { UNIMEAN = 2, LAPLACEBELTRAMI = 3, GAUSS = 4 };

    // Boolean for the viewer
    bool wireframe = false;
    /** NEW: surface display boolean **/
    bool displaySurface = true;
    /** NEW! Temperature boolean**/
    bool temperature = false;
    /** End NEW */

    // Grab variables
    bool m_isGrabbing = false;
    std::vector<Index> m_grabSelection;
    Eigen::Matrix4f m_currentMVP; // Current Model-View-Projection
    float m_grabRadius = INITIAL_GRAB_RADIUS;
    std::vector<Vector3f> m_grabVertPos;
    Vector3f m_grabOrigPos;
    Vector4f m_grabOrigProj;
    Vector2i m_lastMousePosition;
    void (*m_grab_callback)(const std::vector<Index>&,const std::vector<Vector3f>&) = nullptr;
    void (*m_grab_release_callback)() = nullptr;

    // Selection variables
    bool m_isSelecting = false;
    Vector2i m_selectionStart, m_selectionEnd;
    std::vector<Index> m_selectedVertices;

    // Floor height and points
    MatrixXf m_floorPoints;
    float m_floorHeight = 0;
    bool m_floorHeightChanged = false;
    bool m_showFloor = false;
    int m_numFloorFaces = 0;

    PopupButton *popupCurvature;
    FloatBox<float>* coefTextBox;
    IntBox<int>* iterationTextBox;

    // Mesh informations
    int n_vertices = 0;
    int n_faces = 0;
    int n_edges = 0;

    Scalar c_min_value;
    Scalar c_max_value;

	// Temporary storage for updated vertex positions or normals that have to be uploaded to the shader
	// in the next drawContents() call.
	// Can be set using updateShaderVertices() / updateShaderNormals();
	MatrixXf m_updated_shader_verts;
    MatrixXf m_updated_shader_normals;
    MatrixXf m_updated_vertex_selections;

    /** **/
    MatrixXf m_updated_vertex_sources;

    Eigen::Matrix<int, 1, -1> m_current_vertex_status;

    /** NEW! Vertex colors for temperature! **/
    Eigen::MatrixXf m_color_temperature_;

    // Flags that will be set to true when new vertex positions or normals have been povided via updateShaderVertices
	// which need to be re-uploaded at the beginning of the next drawContents() call
	bool m_reupload_vertices = false;
    bool m_reupload_normals = false;
    bool m_reupload_vertex_selections = false;

    /** NEW: Color boolean for reupload purposes **/
    bool m_reupload_vertex_colors = false;

	// Callback function to be called before each draw
	bool (*m_pre_draw_callback)(Viewer*) = nullptr;

	// Callback function to be called after loading a mesh
	bool (*m_mesh_load_callback)(Viewer*) = nullptr;

    // Callback function to be called before loading a mesh, 
    // but after the load button has been pushed
    bool (*m_pre_mesh_load_callback)(Viewer*) = nullptr;


};
